#!/usr/bin/env bash
# Trim adapters with Trim Galore -> crop both reads to fixed length (Cutadapt 3.5 compatible) -> FastQC
# Expects files like name_R1_001.fastq.gz and name_R2_001.fastq.gz

set -euo pipefail

INDIR="."
OUTTG="01_trim_galore"
OUTCROP="02_crop50"
OUTQC="03_fastqc"
TMPDIR="${OUTCROP}/tmp"

THREADS=${THREADS:-8}
CROP_LEN=50       # final per-read length
MINLEN=25         # drop pairs if any read < MINLEN after cropping
NEXTERA=0         # -x to enable Trim Galore --nextera
FORCE=0           # -f to overwrite existing cropped outputs

usage() {
  cat <<EOF
Usage: $(basename "$0") [-i INPUT_DIR] [-t THREADS] [-l CROP_LEN] [-m MINLEN] [-x] [-f]
  -i  Input directory with FASTQs (default: .)
  -t  Threads/cores to use (default: ${THREADS})
  -l  Crop length in bp (default: ${CROP_LEN})
  -m  Minimum length to keep (pair dropped if any read < m) (default: ${MINLEN})
  -x  Use Trim Galore --nextera (default: off)
  -f  Force overwrite existing cropped outputs (default: off)
  -h  Show this help
EOF
  exit 1
}

while getopts ":i:t:l:m:xfh" opt; do
  case $opt in
    i) INDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    l) CROP_LEN="$OPTARG" ;;
    m) MINLEN="$OPTARG" ;;
    x) NEXTERA=1 ;;
    f) FORCE=1 ;;
    h|*) usage ;;
  esac
done

# Dependencies
for x in trim_galore cutadapt fastqc; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: '$x' not found in PATH."; exit 1; }
done

mkdir -p "$OUTTG/logs" "$OUTCROP/logs" "$OUTQC" "$TMPDIR"

shopt -s nullglob
R1S=("${INDIR}"/*_R1_001.fastq.gz)
shopt -u nullglob
[ ${#R1S[@]} -gt 0 ] || { echo "No R1 files found in '$INDIR' matching *_R1_001.fastq.gz"; exit 1; }

echo "Trim Galore : $(trim_galore --version 2>/dev/null || echo unknown)"
echo "cutadapt    : $(cutadapt --version 2>/dev/null || echo unknown) (two-pass crop for Cutadapt 3.5)"
echo "FastQC      : $(fastqc --version 2>/dev/null || echo unknown)"
echo
echo "Input dir   : $INDIR"
echo "Threads     : $THREADS"
echo "Crop length : $CROP_LEN"
echo "Min length  : $MINLEN"
echo "Nextera     : $([ $NEXTERA -eq 1 ] && echo on || echo off)"
echo

for R1 in "${R1S[@]}"; do
  base=$(basename "$R1")
  sample=${base%_R1_001.fastq.gz}
  R2="${INDIR}/${sample}_R2_001.fastq.gz"
  if [ ! -f "$R2" ]; then
    echo "WARNING: Missing R2 for sample '$sample' - skipping."
    continue
  fi

  TG_R1="${OUTTG}/${sample}_R1_001_val_1.fq.gz"
  TG_R2="${OUTTG}/${sample}_R2_001_val_2.fq.gz"

  OUT1="${OUTCROP}/${sample}_R1_001.trim${CROP_LEN}.fastq.gz"
  OUT2="${OUTCROP}/${sample}_R2_001.trim${CROP_LEN}.fastq.gz"
  LOG_TG="${OUTTG}/logs/${sample}.trim_galore.log"
  LOG_CROP1="${OUTCROP}/logs/${sample}.crop${CROP_LEN}.pass1.log"
  LOG_CROP2="${OUTCROP}/logs/${sample}.crop${CROP_LEN}.pass2.log"

  # Skip if outputs already exist and not forcing
  if [ $FORCE -eq 0 ] && [ -s "$OUT1" ] && [ -s "$OUT2" ]; then
    echo ">> [$sample] Cropped outputs exist and --force not set; skipping."
    continue
  fi

  echo ">> [$sample] Trim Galore (adapter/quality trimming)..."
  TG_OPTS=(--paired --gzip --quality 20 --cores "$THREADS" --output_dir "$OUTTG")
  [ $NEXTERA -eq 1 ] && TG_OPTS+=(--nextera)

  if ! trim_galore "${TG_OPTS[@]}" "$R1" "$R2" >"$LOG_TG" 2>&1; then
    echo "ERROR: Trim Galore failed for $sample; see $LOG_TG"; exit 1
  fi
  [[ -s "$TG_R1" && -s "$TG_R2" ]] || { echo "ERROR: Expected Trim Galore outputs not found for $sample"; exit 1; }

  # Temporary outputs after pass 1
  TMP1_R1="${TMPDIR}/${sample}.pass1_R1.fq.gz"  # R1 cropped
  TMP1_R2="${TMPDIR}/${sample}.pass1_R2.fq.gz"  # R2 passthrough

  echo ">> [$sample] Crop pass 1 (crop R1 to ${CROP_LEN}bp; keep pairs synced)..."
  if ! cutadapt \
      -j "$THREADS" \
      -l "$CROP_LEN" \
      --minimum-length "$MINLEN" \
      --pair-filter=any \
      -o "$TMP1_R1" -p "$TMP1_R2" \
      "$TG_R1" "$TG_R2" >"$LOG_CROP1" 2>&1; then
    echo "ERROR: Cutadapt pass 1 failed for $sample; see $LOG_CROP1"; exit 1
  fi

  # Pass 2: swap inputs so original R2 becomes "first read" to crop it
  echo ">> [$sample] Crop pass 2 (crop R2 to ${CROP_LEN}bp by swapping inputs)..."
  if ! cutadapt \
      -j "$THREADS" \
      -l "$CROP_LEN" \
      --minimum-length "$MINLEN" \
      --pair-filter=any \
      -o "$OUT2" -p "$OUT1" \
      "$TMP1_R2" "$TMP1_R1" >"$LOG_CROP2" 2>&1; then
    echo "ERROR: Cutadapt pass 2 failed for $sample; see $LOG_CROP2"; exit 1
  fi

  # Clean up temps
  rm -f "$TMP1_R1" "$TMP1_R2"

  echo ">> [$sample] FastQC on cropped reads..."
  fastqc -t "$THREADS" -o "$OUTQC" "$OUT1" "$OUT2" >/dev/null 2>&1 || {
    echo "ERROR: FastQC failed for $sample"; exit 1; }

  echo ">> [$sample] done."
done

echo
echo "All done:"
echo " - Trimmed (Trim Galore): $OUTTG"
echo " - Cropped to ${CROP_LEN} bp: $OUTCROP"
echo " - FastQC HTML reports: $OUTQC"
