"""Per-day gene count matrices (genes x cells) and joint Leiden labels for sc-SHC."""
import sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.io as sio
J = sys.argv[1]
g = ad.read_h5ad(f"{J}/genes.processed.h5ad")
os.makedirs(f"{J}/scshc", exist_ok=True)
for s in ["day0", "day4"]:
    a = g[(g.obs["sample"] == s).values]
    C = a.layers["counts"].tocsc()
    keep = np.asarray((C > 0).sum(0)).ravel() >= 10
    sio.mmwrite(f"{J}/scshc/{s}_counts.mtx", C[:, keep].T.tocoo())
    pd.Series(a.var_names[keep]).to_csv(f"{J}/scshc/{s}_genes.tsv", index=False, header=False)
    a.obs[["leiden"]].rename_axis("cell").to_csv(f"{J}/scshc/{s}_cells.tsv", sep="\t")
    print(s, a.n_obs, keep.sum())
