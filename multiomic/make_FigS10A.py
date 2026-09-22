"""Fig S10A: GMM BIC across k = 2-10, with a dashed rule at the BIC-best k.
AIC is computed upstream and written to data/gmm_bic_sweep.tsv; it is not drawn.

Output: results/FigS10A/FigS10A.{pdf,png,svg}
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

fs.setup()
DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("FigS10A")

sweep = pd.read_csv(DATA / "gmm_bic_sweep.tsv", sep="\t")
bic_best_k = int(sweep.loc[sweep.BIC.idxmin(), "k"])

fig, ax = plt.subplots(figsize=fs.panel_mm(60, 46))
ax.plot(sweep.k, sweep.BIC, "o-", color="#000000", lw=1.0, ms=3)
ax.axvline(bic_best_k, ls="--", color="#666666", lw=0.8, alpha=0.8)
ax.set_xlabel("k (clusters)", fontsize=8); ax.set_ylabel("BIC", fontsize=8)
ax.tick_params(labelsize=6); fs.prism(ax, lw=1.0)
ax.set_title(f"BIC-best k = {bic_best_k}", fontsize=8)
fig.tight_layout(pad=0.3)
fs.save(fig, OUT / "FigS10A")
print(f"BIC-best k = {bic_best_k}   (sweep k {int(sweep.k.min())}-{int(sweep.k.max())})")
