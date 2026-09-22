#!/bin/bash
# =============================================================================
# Concatenate the H3K27me3 deep-sequencing trim50 FASTQs with the original
# trim50 FASTQs of the same library
# =============================================================================
# gzip files concatenate directly - `cat a.fastq.gz b.fastq.gz > c.fastq.gz`
# produces a valid gzip stream that decoders handle correctly.
#
# Inputs:
#   - NEW_DIR : deep-sequencing fastqs trimmed to 50bp by 00_crop_fastqs.sh
#               (named  5908-0001-Pool-N_S###_L002_R{1,2}_001.trim50.fastq.gz)
#   - OLD_DIR : original fastqs already trimmed to 50bp
#               (named  Naive_H3K27me3_S1_S1_R{1,2}_001.trim50.fastq.gz, etc.)
#
# Output:
#   - OUT_DIR : per-sample concatenated FASTQs named
#               {sample_id}.merged.R{1,2}.fastq.gz
#
# Pool -> sample_id mapping (from the vendor QC summary, project 5908-0001):
#   Pool-1 -> Naive_H3K27me3_R1            (old: Naive_H3K27me3_S1_S1)
#   Pool-2 -> Naive_H3K27me3_R2            (old: Naive_H3K27me3_S2_S2)
#   Pool-3 -> Stim_H3K27me3_R1             (old: Stim_H3K27me3_S1_S3)
#   Pool-4 -> Stim_H3K27me3_R2             (old: Stim_H3K27me3_S2_S4)
#   Pool-5 -> Rested_H3K27me3_R1           (old: Rested_H3K27me3_S1_S5)
#   Pool-6 -> Rested_H3K27me3_R2           (old: Rested_H3K27me3_S2_S6)
#   Pool-7 -> Rechallenged_H3K27me3_R1     (old: Rechallenged_H3K27me3_S1_S7)
#   Pool-8 -> Rechallenged_H3K27me3_R2     (old: Rechallenged_H3K27me3_S2_S8)
# =============================================================================

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

NEW_DIR="${1:-${FASTQ_DIR}/deepseq_trim50/02_crop50}"
OLD_DIR="${2:-${FASTQ_DIR}/GB-249}"
OUT_DIR="${3:-${FASTQ_DIR}/merged}"

mkdir -p "${OUT_DIR}"

# columns: pool_num   new_prefix                      old_prefix                   sample_id
MAP=(
  "1  5908-0001-Pool-1_S112_L002    Naive_H3K27me3_S1_S1            Naive_H3K27me3_R1"
  "2  5908-0001-Pool-2_S113_L002    Naive_H3K27me3_S2_S2            Naive_H3K27me3_R2"
  "3  5908-0001-Pool-3_S114_L002    Stim_H3K27me3_S1_S3             Stim_H3K27me3_R1"
  "4  5908-0001-Pool-4_S115_L002    Stim_H3K27me3_S2_S4             Stim_H3K27me3_R2"
  "5  5908-0001-Pool-5_S116_L002    Rested_H3K27me3_S1_S5           Rested_H3K27me3_R1"
  "6  5908-0001-Pool-6_S117_L002    Rested_H3K27me3_S2_S6           Rested_H3K27me3_R2"
  "7  5908-0001-Pool-7_S118_L002    Rechallenged_H3K27me3_S1_S7     Rechallenged_H3K27me3_R1"
  "8  5908-0001-Pool-8_S119_L002    Rechallenged_H3K27me3_S2_S8     Rechallenged_H3K27me3_R2"
)

for row in "${MAP[@]}"; do
  read -r pool new_pre old_pre sid <<< "${row}"

  for r in R1 R2; do
    new_fq="${NEW_DIR}/${new_pre}_${r}_001.trim50.fastq.gz"
    old_fq="${OLD_DIR}/${old_pre}_${r}_001.trim50.fastq.gz"
    out_fq="${OUT_DIR}/${sid}.merged.${r}.fastq.gz"

    if [[ ! -f "${new_fq}" ]]; then echo "MISSING NEW: ${new_fq}" >&2; exit 1; fi
    if [[ ! -f "${old_fq}" ]]; then echo "MISSING OLD: ${old_fq}" >&2; exit 1; fi

    echo "[MERGE] ${sid} ${r}: ${new_fq} + ${old_fq} -> ${out_fq}"
    cat "${new_fq}" "${old_fq}" > "${out_fq}"

    # Integrity check on the concatenated gzip
    gunzip -t "${out_fq}"
  done
done

echo
echo "Merged FASTQs written to: ${OUT_DIR}"
ls -lh "${OUT_DIR}"
