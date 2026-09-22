#!/bin/bash
# =============================================================================
# Merge replicate BAMs per condition and make the merged CPM bigWigs
# =============================================================================
#   H3K27ac, H3K4me3       - 4 conditions (Naive, Stim, Rested, Rechallenged)
#   H3K36me2, H3K36me3     - 2 conditions (Naive, Stim only)
#   pooled IgG             - single IgG control used for ALL marks (K27/K4 and K36)
# H3K27me3 is merged separately by 04_merge_h3k27me3_deep.sh.
#
#   samtools merge -> samtools sort -> bamCoverage \
#       --normalizeUsing CPM --extendReads --binSize 10 --smoothLength 50
#
# Usage:
#   bash 03_merge_replicates.sh [output_root]
# =============================================================================

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

OUT_ROOT="${1:-./output}"

# Blacklist-filtered BAMs from 02_cutandrun_pipeline.sh (K27ac, K4me3, and the IgG used for ALL marks)
BAM_DIR_GB281="${OUT_ROOT}/GB249/03_filtered"
BAM_SUFFIX_GB281="_hg38_filtered.bam"

# K36me2, K36me3 markdup BAMs from the nf-core/cutandrun run (GB-263 IgG is NOT used)
BAM_DIR_GB263="${K36_BAM_DIR}"
BAM_SUFFIX_GB263=".target.markdup.sorted.bam"

THREADS=16

# Marks x conditions to process here
FULL_CONDITIONS=(Naive Stim Rested Rechallenged)
PARTIAL_CONDITIONS=(Naive Stim)
FULL_MARKS=(H3K27ac H3K4me3)         # all 4 conditions
PARTIAL_MARKS=(K36me2 K36me3)        # Naive + Stim only

mkdir -p "${OUT_ROOT}/merged/merged_bams" "${OUT_ROOT}/merged/bigwigs/cpm"
cd "${OUT_ROOT}/merged"

echo "============================================================"
echo "Replicate merge"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "============================================================"

# =============================================================================
# STEP 0: PREFLIGHT
# =============================================================================
echo ""
echo "--- STEP 0: Preflight checks ---"

all_ok=true

check_bam() {
    local path="$1"
    if [[ -f "$path" ]]; then
        echo "  [OK] $(basename "$path")"
    else
        echo "  [MISSING] $path"
        all_ok=false
    fi
}

# K27ac, K4me3 (4 conds x 2 reps)
for mark in "${FULL_MARKS[@]}"; do
    for cond in "${FULL_CONDITIONS[@]}"; do
        for rep in R1 R2; do
            check_bam "${BAM_DIR_GB281}/${cond}_${mark}_${rep}${BAM_SUFFIX_GB281}"
        done
    done
done

# IgG (2 reps)
for rep in R1 R2; do
    check_bam "${BAM_DIR_GB281}/pooled_IgG_${rep}${BAM_SUFFIX_GB281}"
done

# K36me2, K36me3 (2 conds x 2 reps) - IgG comes from the GB-281 run only
for mark in "${PARTIAL_MARKS[@]}"; do
    for cond in "${PARTIAL_CONDITIONS[@]}"; do
        for rep in R1 R2; do
            check_bam "${BAM_DIR_GB263}/${cond}_${mark}_${rep}${BAM_SUFFIX_GB263}"
        done
    done
done

if [[ "$all_ok" == false ]]; then
    echo ""
    echo "ERROR: one or more required inputs missing. Fix and re-run."
    exit 1
fi

echo "All required inputs present. Proceeding."

# =============================================================================
# STEP 1: MERGE BAMs
# =============================================================================
echo ""
echo "--- STEP 1: Merge BAMs ---"

# Helper: merge R1+R2, sort, index. $1 = output BAM, $2 = R1 BAM, $3 = R2 BAM
merge_pair() {
    local out="$1"
    local r1="$2"
    local r2="$3"
    if [[ -f "$out" ]]; then
        echo "  $out exists, skipping."
        return
    fi
    local unsorted="${out%.bam}_unsorted.bam"
    echo "  Merging $(basename "$r1") + $(basename "$r2") -> $(basename "$out")"
    samtools merge -f -@ "${THREADS}" "$unsorted" "$r1" "$r2"
    samtools sort -@ "${THREADS}" -o "$out" "$unsorted"
    samtools index "$out"
    rm "$unsorted"
    local n=$(samtools view -c "$out")
    echo "    merged reads: $n"
}

for mark in "${FULL_MARKS[@]}"; do
    for cond in "${FULL_CONDITIONS[@]}"; do
        merge_pair \
            "merged_bams/${cond}_${mark}_merged.bam" \
            "${BAM_DIR_GB281}/${cond}_${mark}_R1${BAM_SUFFIX_GB281}" \
            "${BAM_DIR_GB281}/${cond}_${mark}_R2${BAM_SUFFIX_GB281}"
    done
done

# The single IgG control used for ALL marks
merge_pair \
    "merged_bams/pooled_IgG_GB281_merged.bam" \
    "${BAM_DIR_GB281}/pooled_IgG_R1${BAM_SUFFIX_GB281}" \
    "${BAM_DIR_GB281}/pooled_IgG_R2${BAM_SUFFIX_GB281}"

for mark in "${PARTIAL_MARKS[@]}"; do
    for cond in "${PARTIAL_CONDITIONS[@]}"; do
        merge_pair \
            "merged_bams/${cond}_${mark}_merged.bam" \
            "${BAM_DIR_GB263}/${cond}_${mark}_R1${BAM_SUFFIX_GB263}" \
            "${BAM_DIR_GB263}/${cond}_${mark}_R2${BAM_SUFFIX_GB263}"
    done
done

# =============================================================================
# STEP 2: bamCoverage -> CPM bigWigs
# =============================================================================
echo ""
echo "--- STEP 2: Generate merged CPM bigWigs ---"

# Helper: $1 = BAM, $2 = output bigWig, $3 = log
make_bigwig() {
    local bam="$1"
    local out="$2"
    local log="$3"
    if [[ -f "$out" ]]; then
        echo "  $out exists, skipping."
        return
    fi
    echo "  bamCoverage -> $(basename "$out")"
    bamCoverage \
        -b "$bam" \
        --normalizeUsing CPM \
        --extendReads \
        --binSize 10 \
        --smoothLength 50 \
        -p "${THREADS}" \
        -o "$out" \
        2> "$log"
    echo "    done."
}

for mark in "${FULL_MARKS[@]}"; do
    for cond in "${FULL_CONDITIONS[@]}"; do
        make_bigwig \
            "merged_bams/${cond}_${mark}_merged.bam" \
            "bigwigs/cpm/${cond}_${mark}_merged_cpm.bw" \
            "bigwigs/cpm/${cond}_${mark}_merged_bamcoverage_cpm.log"
    done
done

make_bigwig \
    "merged_bams/pooled_IgG_GB281_merged.bam" \
    "bigwigs/cpm/pooled_IgG_GB281_merged_cpm.bw" \
    "bigwigs/cpm/pooled_IgG_GB281_merged_bamcoverage_cpm.log"

for mark in "${PARTIAL_MARKS[@]}"; do
    for cond in "${PARTIAL_CONDITIONS[@]}"; do
        make_bigwig \
            "merged_bams/${cond}_${mark}_merged.bam" \
            "bigwigs/cpm/${cond}_${mark}_merged_cpm.bw" \
            "bigwigs/cpm/${cond}_${mark}_merged_bamcoverage_cpm.log"
    done
done

echo ""
echo "============================================================"
echo "Merge complete: $(date)"
echo "============================================================"
echo ""
echo "Merged BAMs:"
ls -lh merged_bams/*.bam 2>/dev/null | awk '{print "  "$NF, $5}' || echo "  (none)"

echo ""
echo "Merged CPM bigWigs:"
ls bigwigs/cpm/*.bw 2>/dev/null | while read f; do
    echo "  $(basename "$f")"
done || echo "  (none)"
