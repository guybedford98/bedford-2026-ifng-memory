#!/usr/bin/env bash
set -euo pipefail

# Run FiberHMM for the four conditions of one locus.
#
# Usage:
#   bash processing/02_run_fiberhmm_locus_series.sh IDO1
#   bash processing/02_run_fiberhmm_locus_series.sh CXCL11 --force

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash processing/02_run_fiberhmm_locus_series.sh <locus> [--force]"
  exit 2
fi

LOCUS="$(echo "$1" | tr '[:lower:]' '[:upper:]')"
FORCE="${2:-}"
if [[ -n "$FORCE" && "$FORCE" != "--force" ]]; then
  echo "ERROR: optional second argument must be --force"
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$SCRIPT_DIR/../data"
SAMPLE_SETS="$DATA_DIR/locus_sample_sets.tsv"
CONFIGS="$DATA_DIR/locus_configs.tsv"

FIBERHMM_DIR="$(python -c "import sys; sys.path.insert(0, sys.argv[1]); import config; print(config.path('dafseq', 'fiberhmm_dir'))" "$REPO_ROOT/shared")"

for tool in fiberhmm-call fiberhmm-extract; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: $tool was not found in PATH."
    echo "Activate/install FiberHMM before running this series."
    exit 127
  fi
done

if ! awk -F '\t' -v locus="$LOCUS" 'NR > 1 && $1 == locus {found=1} END {exit !found}' "$CONFIGS"; then
  echo "ERROR: locus '$LOCUS' not found in $CONFIGS"
  exit 1
fi

mapfile -t ROWS < <(
  awk -F '\t' -v locus="$LOCUS" 'NR > 1 && $1 == locus {print $4 "\t" $2 "\t" $3}' "$SAMPLE_SETS" | sort -n
)

if [[ "${#ROWS[@]}" -eq 0 ]]; then
  echo "ERROR: no samples configured for locus '$LOCUS' in $SAMPLE_SETS"
  exit 1
fi

echo "=== FiberHMM locus series: $LOCUS ==="
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r _order condition sample <<< "$row"
  out_bam="$FIBERHMM_DIR/$sample.fiberhmm.recalled.bam"
  tf_bed="$FIBERHMM_DIR/$sample.fiberhmm.recalled_tf.bed"

  echo ""
  echo "[$LOCUS] $condition: $sample"
  if [[ "$FORCE" != "--force" && -s "$out_bam" && -s "$tf_bed" ]]; then
    echo "  Existing FiberHMM outputs found; skipping. Use --force to rerun."
    continue
  fi
  bash "$SCRIPT_DIR/01_run_fiberhmm.sh" "$sample"
done

echo ""
echo "$LOCUS FiberHMM series complete. Outputs in $FIBERHMM_DIR"
