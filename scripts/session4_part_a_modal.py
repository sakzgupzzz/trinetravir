"""Modal wrapper for Session 4 Part A scVI sweep.

Bypasses AWS quota / Lambda trust concerns. Runs on Modal (modal.com) — a
serverless GPU platform with pay-per-second billing and auto-terminate on
function return. No instance management; no risk of forgotten billing.

Manual one-time setup:
  1. pip install modal (or: uv add --dev modal)
  2. modal token new → opens browser, authenticate via GitHub/Google
  3. Done. $30 free credit for new accounts covers Session 4 entirely.

Workflow:
  # One-time: create persistent volume + upload inputs
  modal volume create trinetravir-data
  modal volume put trinetravir-data data/processed/scvi_input_monocyte.h5ad /inputs/
  modal volume put trinetravir-data data/processed/scvi_input_B.h5ad        /inputs/
  modal volume put trinetravir-data data/processed/scvi_input_NK.h5ad       /inputs/
  modal volume put trinetravir-data data/processed/scvi_input_CD4T.h5ad     /inputs/
  modal volume put trinetravir-data data/processed/scvi_input_CD8T.h5ad     /inputs/
  modal volume put trinetravir-data data/processed/phase3_response_vectors_monocyte.parquet /inputs/
  modal volume put trinetravir-data data/processed/phase3_response_vectors_B.parquet        /inputs/
  modal volume put trinetravir-data data/processed/phase3_response_vectors_NK.parquet       /inputs/
  modal volume put trinetravir-data data/processed/phase3_response_vectors_CD4T.parquet     /inputs/
  modal volume put trinetravir-data data/processed/phase3_response_vectors_CD8T.parquet     /inputs/

  # Run sweep (auto-launches A10G, runs ~10-12h, auto-terminates)
  modal run scripts/session4_part_a_modal.py

  # Download results
  mkdir -p results/tables
  modal volume get trinetravir-data /outputs/ results/tables/

Cost (Modal pricing 2026-05-11):
  A10G $1.10/hr × 12h = $13.20  (default; sufficient for 4000-HVG scVI)
  A100 $2.78/hr × 7h  = $19.46  (faster; switch GPU="A100" if A10G unavailable)
  T4   $0.59/hr × 18h = $10.62  (cheapest; slower)

Pay-per-second, auto-terminate. Modal $30 new-account credit covers all of these.
"""

from __future__ import annotations

import modal

app = modal.App("trinetravir-session-4")

# Persistent volume for inputs (h5ads, parquets) + outputs (CSVs, logs).
# Volume must be pre-populated with scvi_input_*.h5ad + phase3_response_vectors_*.parquet
# under /inputs/ before `modal run`.
volume = modal.Volume.from_name("trinetravir-data", create_if_missing=True)

# Modal image: CUDA-enabled debian + scientific Python stack pinned to v1
# methods_versions.yaml. uv pip_install respects PyTorch's CUDA wheel index.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential")
    .pip_install(
        "torch==2.11.0",
        extra_index_url="https://download.pytorch.org/whl/cu121",
    )
    .pip_install(
        "scvi-tools==1.4.2",
        "scanpy==1.11.5",
        "anndata==0.12.13",
        "harmonypy==2.0.0",
        "pyarrow",
        "fastparquet",
        "numpy>=1.26,<2.3",
        "pandas>=2.2",
        "scipy>=1.12",
    )
    .add_local_dir(
        local_path=".",
        remote_path="/repo",
        ignore=[
            "**/.git/**",
            "**/.venv/**",
            "**/__pycache__/**",
            "**/data/raw/**",
            "**/data/processed/**",
            "**/data/reference/**.h5ad",
            "**/results/**",
            "**/notebooks/**",
            "**/.terraform/**",
            "**/infrastructure/lambda_state.json",
        ],
    )
)


@app.function(
    image=image,
    gpu="A10G",  # 24GB; sufficient for 4000-HVG scVI on 68K-cell bucket
    timeout=14 * 3600,  # 14h ceiling per spec wall-time estimate + buffer
    volumes={"/data": volume},
)
def run_part_a_sweep() -> None:
    """Run scripts/session4_part_a_scvi_sweep.py with mounted volume."""
    import os
    import subprocess
    import sys

    # Symlink volume inputs into the repo's expected data/processed path
    os.makedirs("/repo/data/processed", exist_ok=True)
    inputs_dir = "/data/inputs"
    if not os.path.isdir(inputs_dir):
        raise FileNotFoundError(
            f"{inputs_dir} not found. Pre-populate via "
            "`modal volume put trinetravir-data data/processed/scvi_input_*.h5ad /inputs/`"
        )
    for fname in os.listdir(inputs_dir):
        src = f"{inputs_dir}/{fname}"
        dst = f"/repo/data/processed/{fname}"
        if not os.path.exists(dst):
            os.symlink(src, dst)
    print(f"Mounted {len(os.listdir(inputs_dir))} input files into /repo/data/processed/")

    # Verify GPU
    import torch

    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU mem: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    # Run sweep
    os.chdir("/repo")
    result = subprocess.run(
        [sys.executable, "scripts/session4_part_a_scvi_sweep.py"],
        check=False,
    )
    print(f"Sweep exit code: {result.returncode}")

    # Copy results into volume /outputs/
    outputs_dir = "/data/outputs"
    os.makedirs(outputs_dir, exist_ok=True)
    import shutil

    src_results = "/repo/results/tables"
    if os.path.isdir(src_results):
        for fname in os.listdir(src_results):
            if fname.startswith("session4_"):
                src = f"{src_results}/{fname}"
                dst = f"{outputs_dir}/{fname}"
                shutil.copy(src, dst)
                print(f"Copied {fname} → /data/outputs/")

    # Persist volume so the user can `modal volume get`
    volume.commit()

    if result.returncode != 0:
        raise RuntimeError(f"Sweep failed with exit code {result.returncode}")


@app.local_entrypoint()
def main() -> None:
    """Launch the sweep. Run via: modal run scripts/session4_part_a_modal.py"""
    print("Launching Session 4 Part A scVI sweep on Modal GPU...")
    print("Expected wall-time: 10-12h on A10G; 7-9h on A100.")
    print("Auto-terminate on function return. Pay-per-second billing.")
    print()
    print("Pre-requisite: volume `trinetravir-data` populated with:")
    print("  /inputs/scvi_input_<bucket>.h5ad × 5")
    print("  /inputs/phase3_response_vectors_<bucket>.parquet × 5")
    print()
    run_part_a_sweep.remote()
    print()
    print("Done. Download results with:")
    print("  modal volume get trinetravir-data /outputs/ results/tables/")
