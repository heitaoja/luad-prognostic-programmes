#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${RSCRIPT:-}" ]]; then
  RSCRIPT="$(command -v Rscript || true)"
  if [[ -z "$RSCRIPT" ]]; then
    RSCRIPT="$(find '/c/Program Files/R' -path '*/bin/Rscript.exe' -print 2>/dev/null | sort -V | tail -1 || true)"
  fi
fi

if [[ -z "$RSCRIPT" ]]; then
  echo "Set RSCRIPT to the Rscript executable." >&2
  exit 1
fi

exec "$RSCRIPT" --vanilla scripts/run_pipeline.R "$@"
