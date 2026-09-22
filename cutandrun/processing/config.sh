#!/bin/bash
# Reference and data paths for the CUT&RUN processing stages, read from the
# repository config.yaml through shared/config.py. Source this, do not run it.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

config_path() {
    python3 -c "import sys; sys.path.insert(0, sys.argv[1]); import config as cfg; print(cfg.path(*sys.argv[2:]))" \
        "${REPO_ROOT}/shared" "$@"
}

HG38_INDEX="$(config_path reference bowtie2_index)"
KMETSTAT_INDEX="$(config_path reference kmetstat_index)"
BLACKLIST="$(config_path reference blacklist)"
SEACR="$(config_path tools seacr)"
FASTQ_DIR="$(config_path cutandrun fastq_dir)"
K36_BAM_DIR="$(config_path cutandrun k36_bam_dir)"
