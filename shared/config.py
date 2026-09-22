"""Reads config.yaml from the repository root and resolves the data paths.

    import sys; sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
    import config as cfg
    counts = cfg.path("rnaseq", "counts")
    out = cfg.out_dir("Fig1E")
"""
import os
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(os.environ.get("IFNG_MSC_CONFIG") or REPO_ROOT / "config.yaml")
CONFIG = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
DATA_ROOT = Path(CONFIG["data_root"]).expanduser()


def path(*keys):
    """Resolve a config key path to an absolute Path, relative keys hanging off data_root."""
    node = CONFIG
    for k in keys:
        node = node[k]
    p = Path(str(node)).expanduser()
    return p if p.is_absolute() else DATA_ROOT / p


def out_dir(*parts):
    """Output directory (created if missing) under the configured output_dir."""
    d = Path(CONFIG["output_dir"]).expanduser()
    if not d.is_absolute():
        d = REPO_ROOT / d
    d = d.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d
