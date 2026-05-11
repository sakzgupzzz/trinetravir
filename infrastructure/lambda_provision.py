"""Lambda Cloud GPU instance provisioning for Session 4 scVI sensitivity.

Bypasses AWS quota restriction by using Lambda Cloud instead. Uses the Lambda
Cloud API (https://cloud.lambdalabs.com/api/v1/docs).

Manual one-time setup (~5 min):
  1. Sign up at https://cloud.lambdalabs.com (use stockbeat email or any).
  2. Add payment method (Dashboard → Account → Billing).
  3. Generate API key (Dashboard → API Keys → "Generate API Key"). Save it.
  4. Upload SSH public key (Dashboard → SSH Keys):
       - Either upload existing ~/.ssh/stockbeat-key.pub
       - Or generate a Lambda-specific key, save private to ~/.ssh/lambda-key.pem
     Name the key reference (e.g., "stockbeat-key" or "lambda-key").

Then export the API key:
  export LAMBDA_API_KEY=secret_xxx...

Usage:
  python infrastructure/lambda_provision.py launch
      → Pick cheapest available GPU instance (A100 SXM4 first; A10 fallback).
      → Launches and waits for "active" state.
      → Writes instance_id + ssh_command to infrastructure/lambda_state.json.

  python infrastructure/lambda_provision.py status
      → Shows current instances + their state.

  python infrastructure/lambda_provision.py terminate
      → Terminates the instance from lambda_state.json.
      → CRITICAL: run this after Session 4 closes to stop billing.

  python infrastructure/lambda_provision.py list-types
      → Shows available instance types + pricing per region.

Cost: A100 SXM4 ~$1.29/hr × 7-9h Session 4 wall-time → ~$9-12.
      A10 24GB ~$0.75/hr × 10-12h → ~$8-9.
      A10 is cheapest absolute; A100 is fastest. Pick A100 first; A10 fallback.

CRITICAL: ALWAYS run `terminate` after Session 4 closes. Lambda keeps billing
until instance is explicitly terminated; no auto-stop. State file at
infrastructure/lambda_state.json tracks the instance for cleanup.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from urllib import error, request

API_BASE = "https://cloud.lambdalabs.com/api/v1"
STATE_FILE = Path(__file__).parent / "lambda_state.json"

# Preference order: fastest + cheapest balance for scVI 200K cells.
PREFERRED_TYPES = [
    "gpu_1x_a100_sxm4",  # $1.29/hr, A100 40GB SXM4
    "gpu_1x_a100",  # $1.29/hr, A100 40GB PCIe
    "gpu_1x_a10",  # $0.75/hr, A10 24GB (sufficient for 4000-HVG scVI)
    "gpu_1x_a6000",  # $0.80/hr, A6000 48GB (fallback)
    "gpu_1x_h100_pcie",  # $2.49/hr, H100 (overkill but works)
]


def get_api_key() -> str:
    key = os.environ.get("LAMBDA_API_KEY", "").strip()
    if not key:
        print("ERROR: set LAMBDA_API_KEY environment variable.", file=sys.stderr)
        print("Get it at https://cloud.lambdalabs.com/api-keys", file=sys.stderr)
        sys.exit(1)
    return key


def api(method: str, path: str, body: dict | None = None) -> dict:
    key = get_api_key()
    auth = f"{key}:".encode()
    import base64

    auth_b64 = base64.b64encode(auth).decode()
    headers = {"Authorization": f"Basic {auth_b64}", "Content-Type": "application/json"}
    data = json.dumps(body).encode() if body else None
    req = request.Request(f"{API_BASE}{path}", data=data, method=method, headers=headers)
    try:
        with request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read())
    except error.HTTPError as e:
        body_str = e.read().decode()
        print(f"HTTP {e.code} on {method} {path}: {body_str}", file=sys.stderr)
        raise


def list_types() -> dict:
    """List instance types + availability + regions."""
    out = api("GET", "/instance-types")
    print(f"{'TYPE':35s} {'PRICE/HR':>10s}  REGIONS_WITH_CAPACITY")
    print("-" * 100)
    types = out.get("data", {})
    for type_name, info in sorted(types.items()):
        regions = info.get("regions_with_capacity_available", [])
        price = info.get("instance_type", {}).get("price_cents_per_hour", 0) / 100
        region_names = ",".join(r["name"] for r in regions) if regions else "(none available)"
        marker = "  <-- preferred" if type_name in PREFERRED_TYPES else ""
        print(f"{type_name:35s} ${price:>8.2f}/hr  {region_names}{marker}")
    return out


def pick_available_type() -> tuple[str, str]:
    """From PREFERRED_TYPES, return the first (type, region) pair with capacity."""
    out = api("GET", "/instance-types")
    types = out.get("data", {})
    for preferred in PREFERRED_TYPES:
        info = types.get(preferred)
        if not info:
            continue
        regions = info.get("regions_with_capacity_available", [])
        if regions:
            region = regions[0]["name"]
            price_cents = info.get("instance_type", {}).get("price_cents_per_hour", 0)
            print(f"Selected: {preferred} in {region} at ${price_cents / 100:.2f}/hr")
            return preferred, region
    raise RuntimeError(
        "No preferred GPU types available in any region. "
        "Run `list-types` to see capacity. Try again later."
    )


def get_ssh_keys() -> list[str]:
    """Return list of ssh key names registered with Lambda Cloud."""
    out = api("GET", "/ssh-keys")
    return [k["name"] for k in out.get("data", [])]


def launch() -> None:
    """Launch a GPU instance + save state."""
    if STATE_FILE.exists():
        state = json.loads(STATE_FILE.read_text())
        print(f"ERROR: state file exists at {STATE_FILE}", file=sys.stderr)
        print(f"  Instance: {state.get('instance_id')} ({state.get('type')})", file=sys.stderr)
        print(
            "  Run `terminate` first if this is a stale state, or `status` to check current state.",
            file=sys.stderr,
        )
        sys.exit(1)

    instance_type, region = pick_available_type()

    ssh_keys = get_ssh_keys()
    if not ssh_keys:
        print("ERROR: no SSH keys registered with Lambda Cloud.", file=sys.stderr)
        print("  Upload via https://cloud.lambdalabs.com/ssh-keys", file=sys.stderr)
        sys.exit(1)
    ssh_key = ssh_keys[0]
    print(f"Using SSH key: {ssh_key}")

    body = {
        "region_name": region,
        "instance_type_name": instance_type,
        "ssh_key_names": [ssh_key],
        "quantity": 1,
        "name": "trinetravir-session-4-gpu",
    }
    print(f"Launching {instance_type} in {region}...")
    out = api("POST", "/instance-operations/launch", body)
    instance_ids = out.get("data", {}).get("instance_ids", [])
    if not instance_ids:
        print(f"ERROR: launch returned no instance_id: {out}", file=sys.stderr)
        sys.exit(1)
    instance_id = instance_ids[0]
    print(f"Instance ID: {instance_id}")

    # Poll until active
    print("Waiting for instance to become active...")
    for attempt in range(60):
        time.sleep(10)
        details = api("GET", f"/instances/{instance_id}")
        status = details.get("data", {}).get("status")
        print(f"  [{attempt + 1}/60] status={status}")
        if status == "active":
            ip = details["data"].get("ip")
            ssh_user = "ubuntu"
            ssh_command = f"ssh {ssh_user}@{ip}"
            STATE_FILE.write_text(
                json.dumps(
                    {
                        "instance_id": instance_id,
                        "type": instance_type,
                        "region": region,
                        "ip": ip,
                        "ssh_user": ssh_user,
                        "ssh_command": ssh_command,
                        "ssh_key_name": ssh_key,
                    },
                    indent=2,
                )
            )
            print()
            print(f"Instance active. IP: {ip}")
            print(f"SSH command:        {ssh_command}")
            print(f"State saved to:     {STATE_FILE}")
            print()
            print("Next steps:")
            print("  1. scp data from local:")
            print(
                f"     scp data/processed/scvi_input_*.h5ad data/processed/phase3_response_vectors_*.parquet {ssh_user}@{ip}:~/"
            )
            print("  2. ssh in + bootstrap:")
            print(f"     {ssh_command}")
            print("     curl -LsSf https://astral.sh/uv/install.sh | sh")
            print("     source $HOME/.local/bin/env")
            print("     git clone https://github.com/sakzgupzzz/trinetravir.git")
            print("     cd trinetravir")
            print("     mkdir -p data/processed && mv ~/scvi_input_*.h5ad data/processed/")
            print("     mv ~/phase3_response_vectors_*.parquet data/processed/")
            print("     uv sync")
            print("     uv run python -c 'import torch; print(torch.cuda.get_device_name(0))'")
            print("  3. Run sweep:")
            print(
                "     nohup uv run python scripts/session4_part_a_scvi_sweep.py > data/processed/session4_part_a.log 2>&1 &"
            )
            print("     tail -f data/processed/session4_part_a.log")
            print("  4. AFTER complete: scp results back + terminate:")
            print(
                f"     scp {ssh_user}@{ip}:~/trinetravir/results/tables/session4_*.csv results/tables/"
            )
            print("     python infrastructure/lambda_provision.py terminate")
            return
    print("ERROR: instance did not become active within 10 min", file=sys.stderr)
    sys.exit(1)


def status() -> None:
    """Show current Lambda instances + state file."""
    if STATE_FILE.exists():
        state = json.loads(STATE_FILE.read_text())
        print(f"Local state file: {STATE_FILE}")
        print(json.dumps(state, indent=2))
        print()

    out = api("GET", "/instances")
    instances = out.get("data", [])
    if not instances:
        print("No instances currently running.")
        return
    print(f"{'INSTANCE_ID':20s} {'TYPE':25s} {'STATUS':12s} {'IP':16s} REGION")
    print("-" * 90)
    for inst in instances:
        print(
            f"{inst['id']:20s} {inst['instance_type']['name']:25s} "
            f"{inst['status']:12s} {inst.get('ip', '-'):16s} {inst['region']['name']}"
        )


def terminate() -> None:
    """Terminate the instance from state file."""
    if not STATE_FILE.exists():
        print(f"ERROR: no state file at {STATE_FILE}.", file=sys.stderr)
        print(
            "  Run `status` to see currently running instances and terminate manually:",
            file=sys.stderr,
        )
        print(
            "  curl -u $LAMBDA_API_KEY: -X POST $API_BASE/instance-operations/terminate \\",
            file=sys.stderr,
        )
        print(
            '       -H "Content-Type: application/json" -d \'{"instance_ids": ["<id>"]}\'',
            file=sys.stderr,
        )
        sys.exit(1)
    state = json.loads(STATE_FILE.read_text())
    instance_id = state["instance_id"]
    print(f"Terminating {instance_id} ({state.get('type')})...")
    out = api("POST", "/instance-operations/terminate", {"instance_ids": [instance_id]})
    print(f"Terminate response: {json.dumps(out, indent=2)}")
    STATE_FILE.unlink()
    print("State file removed. Billing stopped.")


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    cmd = sys.argv[1]
    if cmd == "launch":
        launch()
    elif cmd == "status":
        status()
    elif cmd == "terminate":
        terminate()
    elif cmd == "list-types":
        list_types()
    else:
        print(f"Unknown command: {cmd}", file=sys.stderr)
        print("Valid: launch | status | terminate | list-types", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
