#!/bin/bash
# Runs the CUT&RUN processing stages in order.
#   bash run_all.sh [output_root]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh"

OUT_ROOT="${1:-${here}/output}"

# H3K27me3 deep-sequencing arm: crop the new reads to 50 bp to match the original
# libraries, then concatenate the two sets per library.
mkdir -p "${FASTQ_DIR}/deepseq_trim50"
( cd "${FASTQ_DIR}/deepseq_trim50" && bash "${here}/00_crop_fastqs.sh" -i "${FASTQ_DIR}/deepseq_raw" -t 8 )
bash "${here}/01_merge_resequenced_fastqs.sh"

for sheet in GB249 GB263 GB301; do
    bash "${here}/02_cutandrun_pipeline.sh" "${here}/data/samples_${sheet}.tsv" "${OUT_ROOT}/${sheet}"
done

bash "${here}/03_merge_replicates.sh" "${OUT_ROOT}"
bash "${here}/04_merge_h3k27me3_deep.sh" "${OUT_ROOT}/GB301"
