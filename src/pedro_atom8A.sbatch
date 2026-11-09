#!/bin/bash
#SBATCH --job-name=pedro_atom8A
##SBATCH --account=<your_slurm_account>
##SBATCH --partition=<your_slurm_partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --array=0-69%4
#SBATCH --output=%x_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec .venv/bin/python -u src/recount_pedro_atom8A.py --index "$SLURM_ARRAY_TASK_ID"
