#!/usr/bin/env bash
set -euo pipefail

# Run FiberHMM for one sample from data/sample_manifest.tsv.
#
# Usage:
#   bash processing/01_run_fiberhmm.sh GB-270_IDO1_Naive_merged
#
# FiberHMM must already be installed in the active environment; this script does
# not install packages.

if [[ $# -ne 1 ]]; then
  echo "Usage: bash processing/01_run_fiberhmm.sh <sample>"
  exit 2
fi

SAMPLE="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MANIFEST="$SCRIPT_DIR/../data/sample_manifest.tsv"

cfg() {
  python -c "import sys; sys.path.insert(0, sys.argv[1]); import config; print(config.path(*sys.argv[2:]))" \
    "$REPO_ROOT/shared" "$@"
}

BAM_DIR="$(cfg dafseq bam_dir)"
FIBERHMM_DIR="$(cfg dafseq fiberhmm_dir)"
REFERENCE="$(cfg reference hg38_fasta)"

if ! command -v fiberhmm-call >/dev/null 2>&1; then
  echo "ERROR: fiberhmm-call was not found in PATH."
  echo "Activate/install FiberHMM v2.10.0+ before running."
  exit 127
fi

ROW="$(awk -F '\t' -v sample="$SAMPLE" 'NR > 1 && $1 == sample {print; exit}' "$MANIFEST")"
if [[ -z "$ROW" ]]; then
  echo "ERROR: sample '$SAMPLE' not found in $MANIFEST"
  exit 1
fi

IFS=$'\t' read -r SAMPLE_NAME LOCUS REGION DECORATED_BAM MAPPED_BAM NOTES <<< "$ROW"

INPUT_BAM="$BAM_DIR/$MAPPED_BAM"
OUT_BAM="$FIBERHMM_DIR/$SAMPLE_NAME.fiberhmm.recalled.bam"
CHROM="${REGION%%:*}"

if [[ ! -s "$INPUT_BAM" ]]; then
  echo "ERROR: input BAM not found or empty: $INPUT_BAM"
  exit 1
fi

mkdir -p "$FIBERHMM_DIR"

echo "Sample: $SAMPLE_NAME"
echo "Locus : $LOCUS"
echo "Region: $REGION"
echo "Input : $INPUT_BAM"
echo "Output: $OUT_BAM"

# NUMBA_DISABLE_JIT avoids a Numba import-cache bug in the installed FiberHMM wheel.
NUMBA_DISABLE_JIT=1 fiberhmm-call \
  -c 4 \
  --io-threads 4 \
  --enzyme ddda \
  --reference "$REFERENCE" \
  --region-parallel \
  --chroms "$CHROM" \
  -i "$INPUT_BAM" \
  -o "$OUT_BAM"

NUMBA_DISABLE_JIT=1 fiberhmm-extract \
  -i "$OUT_BAM" \
  -o "$FIBERHMM_DIR" \
  -c 4 \
  --footprint \
  --msp \
  --tf \
  --keep-bed

echo "Done: $FIBERHMM_DIR"
