#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- run the pipeline without `make`
# -----------------------------------------------------------------------------
# Usage:
#   ./run_all.sh              # run every implemented step in order
#   ./run_all.sh 02 03        # run only the steps whose name starts with 02/03
#   RSCRIPT=/usr/bin/Rscript ./run_all.sh
#
# Works in Git Bash / WSL / Linux. On Windows, if Rscript is not on PATH the
# script looks for the newest R installation under "C:/Program Files/R".
# =============================================================================
set -euo pipefail

STEPS=(
  00_setup
  01_load_data
  02_score
  03_tcga_primary
  03b_mutual_adjustment
  04_external_meta
  05_incremental
  06_negative_control
  07_figures
  08_tables
)

# ---- locate Rscript --------------------------------------------------------
if [[ -z "${RSCRIPT:-}" ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="$(command -v Rscript)"
  else
    RSCRIPT="$(ls -d "/c/Program Files/R"/R-*/bin/Rscript.exe 2>/dev/null | sort -V | tail -1 || true)"
  fi
fi
if [[ -z "${RSCRIPT:-}" || ! -x "$RSCRIPT" ]]; then
  echo "ERROR: Rscript not found. Set RSCRIPT=/path/to/Rscript" >&2
  exit 1
fi
echo "Using R: $RSCRIPT"

mkdir -p data/interim data/processed results/logs results/objects results/figures results/tables

# ---- select steps ----------------------------------------------------------
if [[ $# -gt 0 ]]; then
  SELECTED=()
  for want in "$@"; do
    for s in "${STEPS[@]}"; do
      [[ "$s" == "$want"* ]] && SELECTED+=("$s")
    done
  done
else
  SELECTED=("${STEPS[@]}")
fi

# ---- run -------------------------------------------------------------------
start=$(date +%s)
for s in "${SELECTED[@]}"; do
  if [[ ! -f "R/$s.R" ]]; then
    echo ""
    echo "SKIP  $s  (not implemented yet -- see docs/analysis_plan.md)"
    continue
  fi
  echo ""
  echo "=============================================================="
  echo "  $s    $(date '+%H:%M:%S')"
  echo "=============================================================="
  "$RSCRIPT" "R/$s.R" 2>&1 | tee "results/logs/$s.log"
done
end=$(date +%s)
echo ""
echo "Finished in $(( (end-start)/60 )) min $(( (end-start)%60 )) s."
