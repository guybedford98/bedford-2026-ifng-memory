"""Statistics for the bar-graph panels, chosen by experimental design and computed in
the same script that draws the bar (never hand-transcribed).

    import plotstats as ps
    res = ps.welch(ctrl, treat)
    fs.sig_bracket(ax, 0, 1, info["ytop"].max(), ps.pstars(res["p"]))
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats as _sps


def pstars(p, four_tier=True):
    """4-tier GraphPad convention: **** p<1e-4, *** p<1e-3, ** p<0.01, * p<0.05, ns.
    Non-finite / None -> 'ns'."""
    if p is None or not np.isfinite(p):
        return "ns"
    if four_tier and p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def welch(a, b):
    """Two groups -> Welch's unpaired t-test (unequal variance)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    t, p = _sps.ttest_ind(a[np.isfinite(a)], b[np.isfinite(b)], equal_var=False)
    return {"test": "Welch's unpaired t-test", "stat": float(t), "p": float(p)}


def oneway_tukey(groups):
    """k groups, all-vs-all -> one-way ANOVA (omnibus) + Tukey HSD pairwise.
    groups: dict{label: array}. Returns omnibus F/p and pairwise dict
    {(labelA,labelB): padj} for every unordered pair (labels sorted as given)."""
    labels = list(groups)
    arrs = [np.asarray(groups[l], float) for l in labels]
    arrs = [a[np.isfinite(a)] for a in arrs]
    F, p = _sps.f_oneway(*arrs)
    # Tukey via statsmodels MultiComparison so pair->padj mapping is explicit.
    from statsmodels.stats.multicomp import MultiComparison
    endog, grp = [], []
    for l, a in zip(labels, arrs):
        endog.extend(a.tolist())
        grp.extend([l] * a.size)
    mc = MultiComparison(np.asarray(endog, float), np.asarray(grp))
    tuk = mc.tukeyhsd()
    uniq = list(mc.groupsunique)
    pairs = list(itertools.combinations(range(len(uniq)), 2))
    pw = {}
    for (i, j), padj, rej in zip(pairs, tuk.pvalues, tuk.reject):
        pw[(uniq[i], uniq[j])] = float(padj)
    return {"test": "one-way ANOVA + Tukey HSD", "F": float(F), "p_omnibus": float(p),
            "pairwise": pw}


def oneway_dunnett(control, others):
    """Guides/conditions vs ONE shared control -> one-way ANOVA omnibus + Dunnett.
    control: array (the reference). others: dict{label: array}.
    Returns {label: p_vs_control} plus omnibus F/p."""
    labels = list(others)
    ctrl = np.asarray(control, float)
    ctrl = ctrl[np.isfinite(ctrl)]
    arrs = [np.asarray(others[l], float) for l in labels]
    arrs = [a[np.isfinite(a)] for a in arrs]
    F, p = _sps.f_oneway(ctrl, *arrs)
    # rng pinned so the QMC p-value integration is reproducible across runs.
    res = _sps.dunnett(*arrs, control=ctrl, rng=0)
    pw = {l: float(pv) for l, pv in zip(labels, res.pvalue)}
    return {"test": "one-way ANOVA + Dunnett vs control", "F": float(F),
            "p_omnibus": float(p), "pairwise": pw}


def twoway(cells, factorA="A", factorB="B", key_pair=None):
    """2x2 factorial -> two-way ANOVA; report the INTERACTION term p (the
    priming/memory test) plus main effects, and a key pairwise (Welch) if asked.
    cells: dict{(levelA, levelB): array}. key_pair: ((a1,b1),(a2,b2)) cells to
    contrast (e.g. NT+IFN vs targeted+IFN). Returns interaction p, main effects,
    and key pairwise p."""
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    rows = []
    for (la, lb), arr in cells.items():
        for v in np.asarray(arr, float):
            if np.isfinite(v):
                rows.append({"A": str(la), "B": str(lb), "y": float(v)})
    df = pd.DataFrame(rows)
    model = smf.ols("y ~ C(A) * C(B)", data=df).fit()
    aov = sm.stats.anova_lm(model, typ=2)

    def _p(term):
        return float(aov.loc[term, "PR(>F)"]) if term in aov.index else np.nan

    out = {"test": "two-way ANOVA (interaction term)",
           "p_interaction": _p("C(A):C(B)"),
           "p_mainA": _p("C(A)"), "p_mainB": _p("C(B)"),
           "factorA": factorA, "factorB": factorB, "anova": aov}
    if key_pair is not None:
        (c1, c2) = key_pair
        out["key_pair"] = (c1, c2)
        out["p_key"] = welch(cells[c1], cells[c2])["p"]
    return out


def fisher(table, alternative="two-sided"):
    """2x2 contingency -> Fisher exact (two-sided by default), rendered on the
    4-tier asterisk scale. table = [[a,b],[c,d]]."""
    odds, p = _sps.fisher_exact(table, alternative=alternative)
    return {"test": "Fisher exact test (two-sided)", "oddsratio": float(odds), "p": float(p)}


def deseq_padj(padj):
    """RNA-seq TPM bars: asterisks come from the existing DESeq2 (Wald + BH-FDR)
    gene-level padj for the relevant contrast; the bar itself is TPM. Do NOT
    recompute a test on the 3 TPM values."""
    return {"test": "DESeq2 on counts (Wald + BH-FDR)",
            "p": float(padj) if padj is not None and np.isfinite(padj) else np.nan}
