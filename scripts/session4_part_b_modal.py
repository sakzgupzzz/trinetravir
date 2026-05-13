"""Modal wrapper for Session 4 Part B global scVI sweep.

Single scVI training run on full v1 corpus (244,389 cells × 4000 HVG matching
harmony_global_embedding.h5ad input space). 16-config grid.

Pre-requisite uploads to Modal volume trinetravir-data /inputs/:
  scvi_input_global.h5ad  (244,389 × 4000; raw counts; obs[study_id, donor_id,
                           donor_disease_status, coarse])
  harmony_global_embedding.h5ad  (Harmony reference for Δr comparison)
  khatri_mvs_gene_list.csv  (already in /inputs/ from Part A)

Launch:
  modal run scripts/session4_part_b_modal.py

Cost (Modal pricing 2026-05-12):
  L4  $1.05/hr × 6h = $6.30 (default)
  A100 $2.78/hr × 3h = $8.34 (faster)

Auto-terminate on function return. Outputs land in volume /outputs/:
  session4_part_b_global_scvi_per_bucket.csv
  session4_part_b_global_verdict.csv
"""

from __future__ import annotations

import modal

app = modal.App("trinetravir-session-4-part-b")

volume = modal.Volume.from_name("trinetravir-data", create_if_missing=True)

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
    gpu="A100",  # 3× faster than L4 on 244K-cell global scVI; ~5-6h vs ~17h on L4
    timeout=8 * 3600,  # 8h ceiling: ~5-6h A100 estimate + buffer
    volumes={"/data": volume},
)
def run_part_b_global_sweep() -> None:
    """Run scripts/session4_part_b_global_sweep.py with mounted volume."""
    import os
    import subprocess
    import sys

    os.makedirs("/repo/data/processed", exist_ok=True)
    inputs_dir = "/data/inputs"
    if not os.path.isdir(inputs_dir):
        raise FileNotFoundError(f"{inputs_dir} missing")
    for fname in os.listdir(inputs_dir):
        src = f"{inputs_dir}/{fname}"
        dst = f"/repo/data/processed/{fname}"
        if not os.path.exists(dst):
            os.symlink(src, dst)
    print(f"Mounted {len(os.listdir(inputs_dir))} input files into /repo/data/processed/")

    import torch

    print(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU mem: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    os.chdir("/repo")
    result = subprocess.run(
        [sys.executable, "scripts/session4_part_b_global_sweep.py"],
        check=False,
    )
    print(f"Part B sweep exit code: {result.returncode}")

    outputs_dir = "/data/outputs"
    os.makedirs(outputs_dir, exist_ok=True)
    import shutil

    src_results = "/repo/results/tables"
    if os.path.isdir(src_results):
        for fname in os.listdir(src_results):
            if fname.startswith("session4_part_b"):
                src = f"{src_results}/{fname}"
                dst = f"{outputs_dir}/{fname}"
                shutil.copy(src, dst)
                print(f"Copied {fname} → /data/outputs/")

    volume.commit()

    if result.returncode != 0:
        raise RuntimeError(f"Part B sweep failed with exit code {result.returncode}")


@app.local_entrypoint()
def main() -> None:
    print("Launching Session 4 Part B (global scVI) on Modal L4...")
    print("Expected wall-time: ~3-6h on L4; ~2-4h on A100.")
    print("Pay-per-second, auto-terminate on return.")
    print()
    print("Pre-requisite volume contents:")
    print("  /inputs/scvi_input_global.h5ad")
    print("  /inputs/harmony_global_embedding.h5ad")
    print("  /inputs/khatri_mvs_gene_list.csv")
    print()
    run_part_b_global_sweep.remote()
    print()
    print("Done. Download results:")
    print(
        "  modal volume get trinetravir-data /outputs/session4_part_b_global_scvi_per_bucket.csv results/tables/"
    )
    print(
        "  modal volume get trinetravir-data /outputs/session4_part_b_global_verdict.csv results/tables/"
    )
