#!/bin/bash
# =============================================================================
# CUT&RUN pipeline - EpiCypher K-MetStat spike-in normalization
# =============================================================================
# Usage:
#   bash 02_cutandrun_pipeline.sh <samples.tsv> [output_dir]
#
# Sample TSV format (tab-separated, with header):
#   sample_id  mark  condition  replicate  R1_path  R2_path  is_igg
# R1_path/R2_path are relative to cutandrun.fastq_dir in config.yaml.
#
# Outputs:
#   01_trimmed/      - fastp trimmed FASTQs + QC reports
#   02_aligned/      - Bowtie2 BAMs (hg38 and K-MetStat)
#   03_filtered/     - Blacklist-filtered hg38 BAMs
#   04_bigwigs/      - Spike-in normalized + CPM normalized bigWigs
#   05_peaks/        - SEACR (primary) and MACS2 (secondary) peak calls
#   06_qc/           - MultiQC, fragment size plots, spike-in stats table
#
# Requirements:
#   bowtie2, samtools, fastp, fastqc, multiqc, deeptools, bedtools, macs2, SEACR
# =============================================================================

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

OUTDIR="./output"                              # Output directory
THREADS="${THREADS:-4}"                        # low concurrency keeps the WSL2 memory/IO footprint small
SORT_THREADS=2                                 # samtools sort threads (kept low)
SORT_MEM="768M"                                # samtools sort memory PER thread (hard cap; total ~1.5G)
RETRIES="${RETRIES:-2}"                        # retries for transient failures of heavy steps
SPIKE_IN_CONSTANT=10000                        # Normalization constant (C in SF = C/spike-in reads)
KEEP_TRIMMED=false                             # Set to true to keep trimmed FASTQs after alignment

# BigWig parameters
BIN_SIZE=10
SMOOTH_LENGTH=50

# Alignments run serially, so each may use all THREADS.
ALIGN_THREADS="${THREADS}"

# Rscript binary directory - passed as arg 6 to the CUTRUNTools SEACR_1.1.sh
RSCRIPT_BIN=$(dirname "$(which Rscript)")

# Parse arguments
SAMPLES_TSV="${1:?Usage: bash 02_cutandrun_pipeline.sh <samples.tsv> [output_dir]}"
[[ -n "${2:-}" ]] && OUTDIR="$2"

# Preflight: required tools must be on PATH (fail fast, not hours into the run)
for tool in bowtie2 samtools fastp fastqc bamCoverage bamPEFragmentSize macs2 multiqc bedtools bc sort awk Rscript; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required tool '$tool' not found in PATH."; exit 1; }
done

# Preflight: reference files must exist
[[ -f "${HG38_INDEX}.1.bt2" || -f "${HG38_INDEX}.1.bt2l" ]] || { echo "ERROR: hg38 bowtie2 index not found at ${HG38_INDEX}.*.bt2"; exit 1; }
[[ -f "${KMETSTAT_INDEX}.1.bt2" || -f "${KMETSTAT_INDEX}.1.bt2l" ]] || { echo "ERROR: K-MetStat bowtie2 index not found at ${KMETSTAT_INDEX}.*.bt2"; exit 1; }
[[ -f "${BLACKLIST}" ]] || { echo "ERROR: blacklist not found at ${BLACKLIST}"; exit 1; }
[[ -f "${SEACR}" ]] || { echo "ERROR: SEACR script not found at ${SEACR}"; exit 1; }

# Create output directory structure
mkdir -p \
    "${OUTDIR}/01_trimmed" \
    "${OUTDIR}/02_aligned" \
    "${OUTDIR}/03_filtered" \
    "${OUTDIR}/04_bigwigs/spikein" \
    "${OUTDIR}/04_bigwigs/cpm" \
    "${OUTDIR}/05_peaks/seacr" \
    "${OUTDIR}/05_peaks/macs2" \
    "${OUTDIR}/06_qc/fastqc" \
    "${OUTDIR}/06_qc/fragment_sizes" \
    "${OUTDIR}/06_qc/kmetstat_barcodes"

# Initialize spike-in stats table (only write the header if absent, so a resumed
# run keeps the rows already computed for completed samples)
STATS_FILE="${OUTDIR}/06_qc/spike_in_stats.tsv"
[[ -f "${STATS_FILE}" ]] || echo -e "sample_id\tmark\tcondition\treplicate\tis_igg\thg38_reads\ttotal_kmetstat_reads\ttarget_kmetstat_reads\tpct_target_spikein\tscale_factor" > "${STATS_FILE}"

# Initialize per-barcode summary (long format: one row per sample x barcode)
BARCODE_SUMMARY="${OUTDIR}/06_qc/kmetstat_barcode_summary.tsv"
[[ -f "${BARCODE_SUMMARY}" ]] || echo -e "sample_id\tmark\tcondition\treplicate\tis_igg\tbarcode\tbarcode_mark\treads" > "${BARCODE_SUMMARY}"

# =============================================================================
# FUNCTIONS
# =============================================================================

# Retry a command up to RETRIES extra times, so one transient failure does not
# halt a long run. A genuine, repeated failure still propagates (set -e aborts).
retry() {
    local n=0
    until "$@"; do
        n=$((n + 1))
        if [[ ${n} -gt ${RETRIES} ]]; then
            echo "ERROR: command failed after $((RETRIES + 1)) attempts: $*" >&2
            return 1
        fi
        echo "  [RETRY ${n}/${RETRIES}] after failure: $* - syncing, waiting 30s" >&2
        sync || true
        sleep 30
    done
}

run_trim() {
    local sample_id="$1" r1="$2" r2="$3"
    local out_r1="${OUTDIR}/01_trimmed/${sample_id}_R1_trimmed.fastq.gz"
    local out_r2="${OUTDIR}/01_trimmed/${sample_id}_R2_trimmed.fastq.gz"

    echo "[TRIM] ${sample_id}"
    fastp \
        -i "${r1}" -I "${r2}" \
        -o "${out_r1}" -O "${out_r2}" \
        --detect_adapter_for_pe \
        -q 20 -l 25 \
        --cut_right --cut_mean_quality 15 \
        -j "${OUTDIR}/01_trimmed/${sample_id}_fastp.json" \
        -h "${OUTDIR}/01_trimmed/${sample_id}_fastp.html" \
        -w "${THREADS}" \
        2> "${OUTDIR}/01_trimmed/${sample_id}_fastp.log"

    fastqc "${out_r1}" "${out_r2}" \
        -o "${OUTDIR}/06_qc/fastqc" \
        -t 2 \
        2>> "${OUTDIR}/01_trimmed/${sample_id}_fastp.log"
}

align_hg38() {
    local sample_id="$1"
    local r1="${OUTDIR}/01_trimmed/${sample_id}_R1_trimmed.fastq.gz"
    local r2="${OUTDIR}/01_trimmed/${sample_id}_R2_trimmed.fastq.gz"
    local bam_out="${OUTDIR}/02_aligned/${sample_id}_hg38.bam"

    echo "[ALIGN hg38] ${sample_id}"
    bowtie2 \
        --local --very-sensitive-local \
        --no-unal --no-mixed --no-discordant \
        --phred33 -I 10 -X 700 \
        -p "${ALIGN_THREADS}" \
        -x "${HG38_INDEX}" \
        -1 "${r1}" -2 "${r2}" \
        2> "${OUTDIR}/02_aligned/${sample_id}_hg38_bowtie2.log" \
        | samtools view -bS -@ "${SORT_THREADS}" \
        | samtools sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -T "${OUTDIR}/02_aligned/${sample_id}_hg38.sorttmp" -o "${bam_out}" -

    samtools index "${bam_out}"
    samtools quickcheck "${bam_out}" \
        || { echo "ERROR: ${sample_id} hg38 BAM failed quickcheck (alignment/sort error)" >&2; return 1; }
}

align_kmetstat() {
    local sample_id="$1"
    local r1="${OUTDIR}/01_trimmed/${sample_id}_R1_trimmed.fastq.gz"
    local bam_out="${OUTDIR}/02_aligned/${sample_id}_kmetstat.bam"

    echo "[ALIGN K-MetStat] ${sample_id}"
    # Single-end R1 only: K-MetStat FASTA has 22bp barcodes; paired-end mode is
    # incompatible with a reference this short. -L 10 overrides the default seed
    # length of 20 (--very-sensitive-local), which is too long for a 22bp reference.
    # --score-min C,22,0 accepts any alignment scoring >= 22 (a ~11bp perfect match),
    # allowing bowtie2 to find the barcode within the 50bp read via soft-clipping.
    #
    # bowtie2 pipes straight into `samtools sort` (which parses SAM natively): the
    # intermediate `samtools view -bS` stage segfaults on the deeper K-MetStat SAM
    # stream, leaving a truncated, unindexed BAM and no scale factor.
    bowtie2 \
        --local --very-sensitive-local \
        -L 10 --score-min C,22,0 \
        --phred33 \
        -p "${ALIGN_THREADS}" \
        -x "${KMETSTAT_INDEX}" \
        -U "${r1}" \
        2> "${OUTDIR}/02_aligned/${sample_id}_kmetstat_bowtie2.log" \
        | samtools sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -T "${OUTDIR}/02_aligned/${sample_id}_kmetstat.sorttmp" -o "${bam_out}" -

    samtools index "${bam_out}"
    samtools quickcheck "${bam_out}" \
        || { echo "ERROR: ${sample_id} K-MetStat BAM failed quickcheck (alignment/sort error)" >&2; return 1; }
}

filter_blacklist() {
    local sample_id="$1"
    local bam_in="${OUTDIR}/02_aligned/${sample_id}_hg38.bam"
    local bam_out="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"

    echo "[FILTER blacklist] ${sample_id}"
    bedtools intersect \
        -v -abam "${bam_in}" \
        -b "${BLACKLIST}" \
        > "${bam_out}"

    samtools index "${bam_out}"
}

save_barcode_counts() {
    local sample_id="$1" mark="$2" condition="$3" replicate="$4" is_igg="$5"
    local bam_kmetstat="${OUTDIR}/02_aligned/${sample_id}_kmetstat.bam"

    echo "[BARCODE COUNTS] ${sample_id}" >&2

    # Per-sample file: one row per barcode
    local outfile="${OUTDIR}/06_qc/kmetstat_barcodes/${sample_id}_barcode_counts.tsv"
    echo -e "barcode\tbarcode_mark\treads" > "${outfile}"
    samtools idxstats "${bam_kmetstat}" | awk '!/^\*/ && NF>=3 {
        split($1, parts, "_")
        barcode_mark = parts[2]
        print $1 "\t" barcode_mark "\t" $3
    }' >> "${outfile}"

    # Append to master summary
    samtools idxstats "${bam_kmetstat}" | awk \
        -v sid="${sample_id}" -v m="${mark}" \
        -v cond="${condition}" -v rep="${replicate}" -v igg="${is_igg}" \
        '!/^\*/ && NF>=3 {
            split($1, parts, "_")
            barcode_mark = parts[2]
            print sid "\t" m "\t" cond "\t" rep "\t" igg "\t" $1 "\t" barcode_mark "\t" $3
        }' >> "${BARCODE_SUMMARY}"
}

calc_scale_factor() {
    local sample_id="$1" mark="$2" condition="$3" replicate="$4" is_igg="$5"
    local bam_hg38="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"
    local bam_kmetstat="${OUTDIR}/02_aligned/${sample_id}_kmetstat.bam"

    echo "[SCALE FACTOR] ${sample_id}" >&2
    local hg38_reads
    hg38_reads=$(samtools view -c -F 4 "${bam_hg38}")

    local total_kmetstat_reads
    total_kmetstat_reads=$(samtools view -c -F 4 "${bam_kmetstat}")

    # Use mark-specific barcode reads for histone marks (e.g. H3K4me3 -> B038_H3K4me3_A + B041_H3K4me3_B).
    # IgG has no specific PTM target, so use total K-MetStat reads instead.
    local target_reads
    if [[ "${is_igg}" == "true" ]] || [[ "${is_igg}" == "1" ]]; then
        target_reads="${total_kmetstat_reads}"
    else
        target_reads=$(samtools idxstats "${bam_kmetstat}" | awk -v mark="${mark}" \
            '$1 ~ mark {sum += $3} END {print sum+0}')
    fi

    if [[ "${target_reads}" -lt 50 ]]; then
        echo "WARNING: ${sample_id} has <50 target K-MetStat reads (${target_reads} for ${mark}). Spike-in normalization skipped; CPM only." >&2
        echo -e "${sample_id}\t${mark}\t${condition}\t${replicate}\t${is_igg}\t${hg38_reads}\t${total_kmetstat_reads}\t${target_reads}\tNA\tNA" >> "${STATS_FILE}"
        echo "NA"
        return
    fi

    local pct_spikein
    pct_spikein=$(echo "scale=4; (${target_reads} / ${hg38_reads}) * 100" | bc)

    local scale_factor
    scale_factor=$(echo "scale=6; ${SPIKE_IN_CONSTANT} / ${target_reads}" | bc)

    echo -e "${sample_id}\t${mark}\t${condition}\t${replicate}\t${is_igg}\t${hg38_reads}\t${total_kmetstat_reads}\t${target_reads}\t${pct_spikein}\t${scale_factor}" >> "${STATS_FILE}"

    echo "${scale_factor}"
}

make_bigwigs() {
    local sample_id="$1" scale_factor="$2"
    local bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"

    echo "[BIGWIG] ${sample_id}"

    # The two bamCoverage calls run serially and under `timeout`: two parallel
    # bigWig writes to a mounted Windows drive deadlock under WSL, and the timeout
    # turns any residual hang into a clean abort rather than an overnight stall.
    if [[ "${scale_factor}" != "NA" ]]; then
        timeout 60m bamCoverage \
            -b "${bam}" \
            -o "${OUTDIR}/04_bigwigs/spikein/${sample_id}_spikein.bw" \
            --scaleFactor "${scale_factor}" \
            --binSize "${BIN_SIZE}" \
            --smoothLength "${SMOOTH_LENGTH}" \
            --normalizeUsing None \
            --extendReads \
            --ignoreDuplicates \
            -p "${THREADS}" \
            2> "${OUTDIR}/04_bigwigs/spikein/${sample_id}_bamcoverage.log"
    else
        echo "  Skipping spike-in bigWig for ${sample_id} (no spike-in reads)"
    fi

    timeout 60m bamCoverage \
        -b "${bam}" \
        -o "${OUTDIR}/04_bigwigs/cpm/${sample_id}_cpm.bw" \
        --normalizeUsing CPM \
        --binSize "${BIN_SIZE}" \
        --smoothLength "${SMOOTH_LENGTH}" \
        --extendReads \
        --ignoreDuplicates \
        -p "${THREADS}" \
        2> "${OUTDIR}/04_bigwigs/cpm/${sample_id}_bamcoverage_cpm.log"
}

call_peaks_seacr() {
    local sample_id="$1" igg_bedgraph="$2" scale_factor="$3"
    local bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"
    local bedgraph="${OUTDIR}/05_peaks/seacr/${sample_id}.bedgraph"

    echo "[PEAKS SEACR] ${sample_id}"

    # Generate spike-in normalized bedGraph for SEACR input
    local sf_arg
    if [[ "${scale_factor}" != "NA" ]]; then
        sf_arg="${scale_factor}"
    else
        sf_arg="1"  # Fall back to unscaled
    fi

    timeout 90m bamCoverage \
        -b "${bam}" \
        -o "${bedgraph}" \
        --outFileFormat bedgraph \
        --scaleFactor "${sf_arg}" \
        --normalizeUsing None \
        --binSize 1 \
        --extendReads \
        --ignoreDuplicates \
        -p "${THREADS}" \
        2>/dev/null

    # Sort bedGraph (SEACR requires sorted input)
    sort -k1,1 -k2,2n "${bedgraph}" > "${bedgraph}.sorted"
    mv "${bedgraph}.sorted" "${bedgraph}"

    # Run SEACR (stringent and relaxed)
    # R_LIBS_USER is unset to prevent an older personal R library from shadowing
    # the conda env's packages (IRanges version conflict with GenomicRanges).
    R_LIBS_USER="" bash "${SEACR}" \
        "${bedgraph}" "${igg_bedgraph}" \
        non stringent \
        "${OUTDIR}/05_peaks/seacr/${sample_id}_seacr" \
        "${RSCRIPT_BIN}" \
        2> "${OUTDIR}/05_peaks/seacr/${sample_id}_seacr.log"

    R_LIBS_USER="" bash "${SEACR}" \
        "${bedgraph}" "${igg_bedgraph}" \
        non relaxed \
        "${OUTDIR}/05_peaks/seacr/${sample_id}_seacr_relaxed" \
        "${RSCRIPT_BIN}" \
        2>> "${OUTDIR}/05_peaks/seacr/${sample_id}_seacr.log"
}

call_peaks_macs2() {
    local sample_id="$1" igg_bam="$2"
    local bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"

    echo "[PEAKS MACS2] ${sample_id}"
    macs2 callpeak \
        -t "${bam}" \
        -c "${igg_bam}" \
        -f BAMPE -g hs \
        --keep-dup all \
        --nomodel \
        -n "${sample_id}" \
        --outdir "${OUTDIR}/05_peaks/macs2/" \
        2> "${OUTDIR}/05_peaks/macs2/${sample_id}_macs2.log"
}

fragment_size_qc() {
    local sample_id="$1"
    local bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"

    # Non-fatal (|| true) so a QC-only hang cannot kill the run
    timeout 30m bamPEFragmentSize \
        -b "${bam}" \
        --histogram "${OUTDIR}/06_qc/fragment_sizes/${sample_id}_fragment_sizes.png" \
        --plotTitle "${sample_id} Fragment Sizes" \
        -p "${THREADS}" \
        2>/dev/null || echo "  WARNING: fragment_size_qc skipped/failed for ${sample_id}"
}

# =============================================================================
# MAIN LOOP
# =============================================================================

echo "============================================="
echo " CUT&RUN pipeline"
echo " Samples: ${SAMPLES_TSV}"
echo " Output:  ${OUTDIR}"
echo "============================================="

# Arrays to collect IgG info for peak calling
declare -A IGG_SCALE_FACTOR  # sample_id -> scale_factor
declare -A IGG_BAM            # sample_id -> filtered BAM path
declare -A IGG_BEDGRAPH       # sample_id -> bedgraph path

# First pass: process all samples, collect IgG info
while IFS=$'\t' read -r sample_id mark condition replicate r1_path r2_path is_igg; do
    [[ "${sample_id}" == "sample_id" ]] && continue  # Skip header
    [[ -z "${sample_id}" ]] && continue

    r1_path="${FASTQ_DIR}/${r1_path}"
    r2_path="${FASTQ_DIR}/${r2_path}"

    echo ""
    echo "--- Processing: ${sample_id} (${mark}, ${condition}, Rep${replicate}) ---"

    # Check input FASTQs exist
    if [[ ! -f "${r1_path}" ]] || [[ ! -f "${r2_path}" ]]; then
        echo "ERROR: FASTQ files not found for ${sample_id}:"
        echo "  R1: ${r1_path}"
        echo "  R2: ${r2_path}"
        echo "Skipping ${sample_id}."
        continue
    fi

    # -------------------------------------------------------------------------
    # Skip a sample that is already fully processed. "Done" means the filtered BAM
    # and CPM bigWig exist, and (for non-IgG) either the spike-in bigWig exists or
    # the stats row recorded NA (CPM-only), so a re-run picks up where an
    # interrupted run left off.
    # -------------------------------------------------------------------------
    cpm_bw="${OUTDIR}/04_bigwigs/cpm/${sample_id}_cpm.bw"
    spikein_bw="${OUTDIR}/04_bigwigs/spikein/${sample_id}_spikein.bw"
    filtered_bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"
    if [[ -s "${filtered_bam}" && -s "${cpm_bw}" ]]; then
        if [[ "${is_igg}" == "true" || "${is_igg}" == "1" || -s "${spikein_bw}" ]] \
           || grep -qP "^${sample_id}\t.*\tNA$" "${STATS_FILE}" 2>/dev/null; then
            echo "--- SKIP (already complete): ${sample_id} ---"
            if [[ "${is_igg}" == "true" || "${is_igg}" == "1" ]]; then
                IGG_SCALE_FACTOR["${sample_id}"]=$(grep -P "^${sample_id}\t" "${STATS_FILE}" | cut -f10)
                IGG_BAM["${sample_id}"]="${filtered_bam}"
                IGG_BEDGRAPH["${sample_id}"]="${OUTDIR}/05_peaks/seacr/${sample_id}.bedgraph"
            fi
            continue
        fi
    fi

    # Not fully done: scrub any partial rows from a prior interrupted run so this
    # sample's stats/barcode entries are written exactly once when reprocessed.
    if grep -qP "^${sample_id}\t" "${STATS_FILE}" 2>/dev/null; then
        grep -vP "^${sample_id}\t" "${STATS_FILE}" > "${STATS_FILE}.tmp" && mv "${STATS_FILE}.tmp" "${STATS_FILE}"
    fi
    if grep -qP "^${sample_id}\t" "${BARCODE_SUMMARY}" 2>/dev/null; then
        grep -vP "^${sample_id}\t" "${BARCODE_SUMMARY}" > "${BARCODE_SUMMARY}.tmp" && mv "${BARCODE_SUMMARY}.tmp" "${BARCODE_SUMMARY}"
    fi

    run_trim "${sample_id}" "${r1_path}" "${r2_path}"
    # Alignments run serially (one aligner + one sort at a time) to keep the
    # instantaneous memory/IO footprint small, each retried on transient failure.
    retry align_hg38 "${sample_id}"
    retry align_kmetstat "${sample_id}"

    # Delete trimmed FASTQs - no longer needed after both alignments
    if [[ "${KEEP_TRIMMED}" == "false" ]]; then
        rm -f "${OUTDIR}/01_trimmed/${sample_id}_R1_trimmed.fastq.gz" \
              "${OUTDIR}/01_trimmed/${sample_id}_R2_trimmed.fastq.gz"
    fi

    filter_blacklist "${sample_id}"

    # Delete unfiltered hg38 BAM - superseded by blacklist-filtered BAM
    rm -f "${OUTDIR}/02_aligned/${sample_id}_hg38.bam" \
          "${OUTDIR}/02_aligned/${sample_id}_hg38.bam.bai"

    save_barcode_counts "${sample_id}" "${mark}" "${condition}" "${replicate}" "${is_igg}"
    scale_factor=$(calc_scale_factor "${sample_id}" "${mark}" "${condition}" "${replicate}" "${is_igg}")

    # Delete K-MetStat BAM - barcode counts and scale factor already recorded
    rm -f "${OUTDIR}/02_aligned/${sample_id}_kmetstat.bam" \
          "${OUTDIR}/02_aligned/${sample_id}_kmetstat.bam.bai"

    make_bigwigs "${sample_id}" "${scale_factor}"
    fragment_size_qc "${sample_id}"

    # Track IgG samples for peak calling
    if [[ "${is_igg}" == "true" ]] || [[ "${is_igg}" == "1" ]]; then
        IGG_SCALE_FACTOR["${sample_id}"]="${scale_factor}"
        IGG_BAM["${sample_id}"]="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"
        IGG_BEDGRAPH["${sample_id}"]="${OUTDIR}/05_peaks/seacr/${sample_id}.bedgraph"

        # Pre-generate IgG bedGraph (needed before peak calling of other samples)
        local_sf="${scale_factor}"
        [[ "${local_sf}" == "NA" ]] && local_sf="1"
        bam="${OUTDIR}/03_filtered/${sample_id}_hg38_filtered.bam"
        bedgraph="${OUTDIR}/05_peaks/seacr/${sample_id}.bedgraph"
        timeout 90m bamCoverage \
            -b "${bam}" \
            -o "${bedgraph}" \
            --outFileFormat bedgraph \
            --scaleFactor "${local_sf}" \
            --normalizeUsing None \
            --binSize 1 \
            --extendReads \
            --ignoreDuplicates \
            -p "${THREADS}" \
            2>/dev/null
        sort -k1,1 -k2,2n "${bedgraph}" > "${bedgraph}.sorted"
        mv "${bedgraph}.sorted" "${bedgraph}"
    fi

    # Flush dirty pages after each sample so the kernel page cache cannot
    # accumulate into a host-starving balloon under WSL2.
    sync

done < "${SAMPLES_TSV}"

# Second pass: call peaks for non-IgG samples
# Uses the first IgG sample found as control (or specify per-sample IgG in TSV)
if [[ "${#IGG_BAM[@]}" -eq 0 ]]; then
    echo "WARNING: No IgG samples found in ${SAMPLES_TSV}. Skipping peak calling."
else
    # Pick first IgG as control
    IGG_CONTROL_ID="${!IGG_BAM[@]}"
    IGG_CONTROL_ID="${IGG_CONTROL_ID%% *}"   # First element
    IGG_CONTROL_BAM="${IGG_BAM[$IGG_CONTROL_ID]}"
    IGG_CONTROL_BEDGRAPH="${IGG_BEDGRAPH[$IGG_CONTROL_ID]}"

    echo ""
    echo "--- Peak Calling (IgG control: ${IGG_CONTROL_ID}) ---"

    while IFS=$'\t' read -r sample_id mark condition replicate r1_path r2_path is_igg; do
        [[ "${sample_id}" == "sample_id" ]] && continue
        [[ -z "${sample_id}" ]] && continue
        if [[ "${is_igg}" == "true" ]] || [[ "${is_igg}" == "1" ]]; then continue; fi  # Skip IgG itself

        # Get previously computed scale factor from stats file (column 10)
        scale_factor=$(grep "^${sample_id}" "${STATS_FILE}" | cut -f10)

        call_peaks_seacr "${sample_id}" "${IGG_CONTROL_BEDGRAPH}" "${scale_factor}"
        call_peaks_macs2 "${sample_id}" "${IGG_CONTROL_BAM}"

        # Delete per-sample bedGraph - large (binSize 1), not needed after peak calling
        rm -f "${OUTDIR}/05_peaks/seacr/${sample_id}.bedgraph"

    done < "${SAMPLES_TSV}"

    # Delete IgG bedGraph - no longer needed after all peak calls complete
    rm -f "${IGG_CONTROL_BEDGRAPH}"
fi

# Final QC: MultiQC
echo ""
echo "--- MultiQC ---"
multiqc "${OUTDIR}" \
    -o "${OUTDIR}/06_qc" \
    -n "multiqc_report" \
    --force \
    2> "${OUTDIR}/06_qc/multiqc.log"

echo ""
echo "============================================="
echo " Pipeline complete"
echo " Spike-in stats: ${STATS_FILE}"
echo " QC report:      ${OUTDIR}/06_qc/multiqc_report.html"
echo "============================================="
