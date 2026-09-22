#!/bin/bash
# =============================================================================
# Merge the deep-sequenced H3K27me3 replicates (R1+R2) per condition
# -> merged bigWigs (CPM + spike-in)
#
# Usage:
#   bash 04_merge_h3k27me3_deep.sh [output_dir_of_the_GB301_run]
# =============================================================================
set -euo pipefail

OUT="${1:-./output/GB301}"
FILT="${OUT}/03_filtered"
STATS="${OUT}/06_qc/spike_in_stats.tsv"
BW_DEST="${OUT}/08_merged/bigwigs"
THREADS=4; SORT_THREADS=2
CONDS=(Naive Stim Rested Rechallenged)

echo "===== H3K27me3 merge - $(date) ====="

command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }
command -v bamCoverage >/dev/null 2>&1 || { echo "ERROR: bamCoverage (deeptools) not found"; exit 1; }
command -v bc >/dev/null 2>&1 || { echo "ERROR: bc not found"; exit 1; }
[[ -f "${STATS}" ]] || { echo "ERROR: ${STATS} not found"; exit 1; }

mkdir -p "${BW_DEST}"

get_target() { grep -P "^$1\t" "${STATS}" | cut -f8; }   # col 8 = target_kmetstat_reads

for cond in "${CONDS[@]}"; do
    r1="${FILT}/${cond}_H3K27me3_R1_hg38_filtered.bam"
    r2="${FILT}/${cond}_H3K27me3_R2_hg38_filtered.bam"
    [[ -f "${r1}" && -f "${r2}" ]] || { echo "ERROR: missing filtered BAM(s) for ${cond}"; exit 1; }
    merged="${BW_DEST}/${cond}_H3K27me3_merged.bam"

    echo ""
    echo "[MERGE] ${cond}  (R1 + R2 -> ${merged##*/})"
    # samtools merge of two coordinate-sorted BAMs yields a sorted BAM (no re-sort needed)
    samtools merge -f -@ "${SORT_THREADS}" "${merged}" "${r1}" "${r2}"
    samtools index "${merged}"
    samtools quickcheck "${merged}" || { echo "ERROR: merged BAM failed quickcheck for ${cond}"; exit 1; }

    # merged spike-in scale factor = 10000 / (target_R1 + target_R2)
    t1=$(get_target "${cond}_H3K27me3_R1"); t2=$(get_target "${cond}_H3K27me3_R2")
    sf=$(echo "scale=6; 10000 / (${t1} + ${t2})" | bc)
    printf "  K-MetStat target reads: %s + %s = %s  ->  merged SF = %s\n" "${t1}" "${t2}" "$((t1 + t2))" "${sf}"

    echo "[BIGWIG cpm]     ${cond}"
    bamCoverage -b "${merged}" -o "${BW_DEST}/${cond}_H3K27me3_merged_cpm.bw" \
        --normalizeUsing CPM --binSize 10 --smoothLength 50 \
        --extendReads --ignoreDuplicates -p "${THREADS}" 2>/dev/null

    echo "[BIGWIG spikein] ${cond}  (SF ${sf})"
    bamCoverage -b "${merged}" -o "${BW_DEST}/${cond}_H3K27me3_merged_spikein.bw" \
        --scaleFactor "${sf}" --normalizeUsing None --binSize 10 --smoothLength 50 \
        --extendReads --ignoreDuplicates -p "${THREADS}" 2>/dev/null

    rm -f "${merged}" "${merged}.bai"   # merged BAM no longer needed (bigWigs made)
    sync
done

echo ""
echo "===== done - $(date) ====="
ls -lh "${BW_DEST}"
