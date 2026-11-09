#!/bin/bash
# Run inside tmux; computation occurs exclusively on Slurm compute nodes.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
out=results/local_interface_geometry/pedro
exec >> "$out/logs/controller.log" 2>> "$out/logs/controller.err"
trap 'rc=$?; if (( rc != 0 )); then echo "Controller stopped with exit $rc; inspect parse/report logs and checkpoints. Resume by rerunning this script in tmux." >&2; fi' EXIT
.venv/bin/python -c 'import json; from pathlib import Path; p=Path("results/local_interface_geometry/pedro"); assert json.loads((p/"validation/status.json").read_text())["status"]=="passed"'
printf 'Controller started: %s\n' "$(date -Is)"
# --wait blocks only this detached tmux controller, never the Codex/SSH foreground.
sbatch --parsable --wait src/pedro_local_geometry.sbatch | tee -a "$out/run_receipt.txt"
sbatch --parsable --wait src/pedro_local_geometry_report.sbatch | tee -a "$out/run_receipt.txt"
printf 'Controller finished: %s\n' "$(date -Is)"
