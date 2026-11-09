#!/bin/bash
#SBATCH --job-name=contact_type_report
##SBATCH --account=<your_slurm_account>
##SBATCH --partition=<your_slurm_partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/contact_type_mpl_$SLURM_JOB_ID
exec .venv/bin/python -u src/report_contact_type_decomposition.py
