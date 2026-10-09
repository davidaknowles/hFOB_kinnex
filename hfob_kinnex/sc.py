"""Single-cell QC, embedding and isoform-usage helpers built on scanpy."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad


def qc_metrics(g):
    """Add per-cell UMI, gene and mitochondrial fraction metrics to a gene-level AnnData."""
    g.var["mt"] = g.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(g, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return g


def mad_outlier(x, nmads=3, side="both"):
    """Flag values more than nmads median absolute deviations from the median."""
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    lo, hi = x < med - nmads * mad, x > med + nmads * mad
    return {"both": lo | hi, "low": lo, "high": hi}[side]


def filter_cells(g, nmads=3, max_mt=None, by="sample"):
    """Per-sample MAD filters on log UMIs, log genes (low side) and mitochondrial percent (high side)."""
    keep = np.ones(g.n_obs, bool)
    for s, idx in g.obs.groupby(by).indices.items():
        o = g.obs.iloc[idx]
        bad = (mad_outlier(np.log1p(o.total_counts.values), nmads, "low")
               | mad_outlier(np.log1p(o.n_genes_by_counts.values), nmads, "low")
               | mad_outlier(o.pct_counts_mt.values, nmads, "high"))
        if max_mt is not None:
            bad |= o.pct_counts_mt.values > max_mt
        keep[idx] = ~bad
    g.obs["qc_pass"] = keep
    return keep


def embed(a, n_hvg=3000, n_pcs=30, resolution=0.5, batch_key=None, scale=True, seed=0):
    """Normalise, log, HVG, PCA, kNN, UMAP and Leiden on a counts AnnData (in place, counts kept in layers)."""
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts", batch_key=batch_key)
    b = a[:, a.var.highly_variable].copy()
    if scale:
        sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=n_pcs, random_state=seed)
    a.obsm["X_pca"] = b.obsm["X_pca"]
    sc.pp.neighbors(a, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2)
    return a


def usage_features(iso, min_gene_cells=0.2, min_iso_frac=0.05, max_genes=3000):
    """Cell x isoform usage matrix (isoform / gene counts), centred per isoform; cells with no gene counts get 0.

    Restricted to genes detected in at least `min_gene_cells` of cells with >=2 isoforms of pooled usage
    >= `min_iso_frac`, keeping the `max_genes` most expressed such genes.
    """
    X = iso.X.tocsc()
    genes = iso.var.gene.values
    tot = np.asarray(X.sum(0)).ravel()
    df = pd.DataFrame({"gene": genes, "tot": tot, "j": np.arange(iso.n_vars)})
    gtot = df.groupby("gene").tot.transform("sum")
    df["frac"] = np.where(gtot > 0, df.tot / gtot, 0)
    df = df[df.frac >= min_iso_frac]
    ok = df.groupby("gene").j.transform("size") >= 2
    df = df[ok]
    gcounts = {}
    rows = []
    for gname, d in df.groupby("gene"):
        gc = np.asarray(X[:, iso.var.gene.values == gname].sum(1)).ravel()
        if (gc > 0).mean() < min_gene_cells:
            continue
        gcounts[gname] = gc
        rows.append((gname, gc.sum()))
    top = [r[0] for r in sorted(rows, key=lambda r: -r[1])[:max_genes]]
    cols, names = [], []
    for gname in top:
        gc = gcounts[gname]
        d = df[df.gene == gname]
        for j, frac in zip(d.j, d.frac):
            x = np.asarray(X[:, j].todense()).ravel()
            u = np.divide(x, gc, out=np.full_like(x, np.nan, dtype=float), where=gc > 0)
            cols.append(np.nan_to_num(u - frac, nan=0.0))
            names.append(iso.var_names[j])
    U = np.column_stack(cols).astype(np.float32)
    return ad.AnnData(X=U, obs=iso.obs.copy(), var=iso.var.loc[names].copy())


def embed_dense(a, n_pcs=30, resolution=0.5, seed=0):
    """PCA/kNN/UMAP/Leiden for an already-transformed dense matrix."""
    sc.tl.pca(a, n_comps=n_pcs, random_state=seed)
    sc.pp.neighbors(a, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2)
    return a
