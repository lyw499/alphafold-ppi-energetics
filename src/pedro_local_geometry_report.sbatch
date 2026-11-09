#!/bin/bash
#SBATCH --job-name=pedro_geometry_report
##SBATCH --account=<your_slurm_account>
##SBATCH --partition=<your_slurm_partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec .venv/bin/python -u src/report_pedro_local_geometry.py
