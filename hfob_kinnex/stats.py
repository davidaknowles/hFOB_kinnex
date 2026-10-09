"""Two-group tests for proportions (isoform / event usage) with cells as observations."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _colsum(M):
    return np.asarray(M.sum(0)).ravel()


def quasibinomial_two_group(Y, N, group, min_cells=20, min_total=50, floor_dispersion=True):
    """Quasi-binomial Wald test of a difference in proportion between two groups, for each column.

    Y : cells x features sparse matrix of successes (e.g. inclusion / isoform counts)
    N : cells x features sparse matrix of trials (e.g. event or gene totals), Y <= N elementwise
    group : boolean array over cells (True = second group)
    The model is logit(p_cg) = a + b * group_c with Var(y) = phi * n p (1-p), equivalent to a quasi-binomial GLM
    with a single binary covariate. phi is estimated by Pearson chi2 over cells with n > 0.
    """
    Y, N = sp.csr_matrix(Y, dtype=np.float64), sp.csr_matrix(N, dtype=np.float64)
    Ninv = N.copy()
    Ninv.data = 1.0 / Ninv.data
    Y2n = Y.multiply(Y).multiply(Ninv).tocsr()
    res = {}
    chi2 = 0.0
    for k, m in (("0", ~group), ("1", group)):
        sy, sn, sy2n = _colsum(Y[m]), _colsum(N[m]), _colsum(Y2n[m])
        ncell = np.asarray((N[m] > 0).sum(0)).ravel()
        p = (sy + 0.5) / (sn + 1.0)
        chi2 = chi2 + (sy2n - 2 * p * sy + p ** 2 * sn) / (p * (1 - p))
        res.update({f"y{k}": sy, f"n{k}": sn, f"cells{k}": ncell, f"p{k}": p})
    nobs = res["cells0"] + res["cells1"]
    phi = chi2 / np.maximum(nobs - 2, 1)
    if floor_dispersion:
        phi = np.maximum(phi, 1.0)
    p0, p1 = res["p0"], res["p1"]
    b = np.log(p1 / (1 - p1)) - np.log(p0 / (1 - p0))
    se = np.sqrt(phi * (1 / (res["n0"] * p0 * (1 - p0)) + 1 / (res["n1"] * p1 * (1 - p1))))
    z = b / se
    df = pd.DataFrame(res)
    df["dPSI"] = df.y1 / np.maximum(df.n1, 1) - df.y0 / np.maximum(df.n0, 1)
    df["logOR"], df["phi"], df["z"] = b, phi, z
    df["pval"] = 2 * stats.norm.sf(np.abs(z))
    testable = (df.cells0 >= min_cells) & (df.cells1 >= min_cells) & (df.n0 >= min_total) & (df.n1 >= min_total)
    df.loc[~testable, ["z", "pval"]] = np.nan
    df["testable"] = testable
    df["padj"] = np.nan
    df.loc[testable, "padj"] = multipletests(df.loc[testable, "pval"], method="fdr_bh")[1]
    return df
