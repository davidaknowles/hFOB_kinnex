"""Cell QC, filtering and embeddings (gene expression, isoform expression, isoform usage) for joint hFOB Kinnex data."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from plotnine import *
from hfob_kinnex.sc import qc_metrics, filter_cells, embed, usage_features, embed_dense

J, REF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
FIG = f"{OUT}/figures"; os.makedirs(FIG, exist_ok=True); os.makedirs(f"{OUT}/tables", exist_ok=True)

g = ad.read_h5ad(f"{J}/genes.h5ad")
iso = ad.read_h5ad(f"{J}/isoforms.h5ad")
qc_metrics(g)
g.obs["n_isoforms"] = np.asarray((iso.X > 0).sum(1)).ravel()
known = iso.var.structural_category.isin(["full-splice_match", "incomplete-splice_match"]).values
g.obs["frac_known_isoform_reads"] = np.asarray(iso.X[:, known].sum(1)).ravel() / g.obs.total_counts
keep = filter_cells(g, nmads=3)
o = g.obs.copy()
o.to_csv(f"{OUT}/tables/cell_qc.tsv", sep="\t")
summ = o.groupby("sample").agg(cells=("qc_pass", "size"), cells_pass=("qc_pass", "sum"),
                               median_umis=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"),
                               median_isoforms=("n_isoforms", "median"), median_pct_mt=("pct_counts_mt", "median"))
summ.to_csv(f"{OUT}/tables/cell_qc_summary.tsv", sep="\t")
print(summ)

m = o.melt(id_vars=["sample", "qc_pass"], value_vars=["total_counts", "n_genes_by_counts", "n_isoforms", "pct_counts_mt"])
p = (ggplot(m, aes("sample", "value", fill="sample")) + geom_violin(scale="width") + geom_boxplot(width=0.1, outlier_size=0)
     + facet_wrap("~variable", scales="free_y", nrow=1) + theme_bw() + theme(legend_position="none") + labs(x="", y=""))
p.save(f"{FIG}/cell_qc_violin.png", width=11, height=3.2, dpi=150)
p = (ggplot(o, aes("total_counts", "n_genes_by_counts", color="qc_pass")) + geom_point(size=0.2, alpha=0.4)
     + facet_wrap("~sample") + scale_x_log10() + scale_y_log10() + theme_bw() + labs(x="UMIs", y="genes"))
p.save(f"{FIG}/cell_qc_scatter.png", width=8, height=3.5, dpi=150)

g = g[keep].copy(); iso = iso[keep].copy()
# drop features with no counts after filtering
g = g[:, np.asarray(g.X.sum(0)).ravel() > 0].copy()
iso = iso[:, np.asarray(iso.X.sum(0)).ravel() > 0].copy()

# gene-level embedding
embed(g, n_hvg=3000, n_pcs=30, resolution=0.5)
cc = pd.read_csv(f"{REF}/cc_genes_seurat2019.tsv", sep="\t")
sc.tl.score_genes_cell_cycle(g, s_genes=[x for x in cc[cc.phase == "S"].gene if x in g.var_names],
                             g2m_genes=[x for x in cc[cc.phase == "G2M"].gene if x in g.var_names])
g.write_h5ad(f"{J}/genes.processed.h5ad")

# isoform-level embedding (log-normalised isoform counts)
iso_e = iso[:, np.asarray((iso.X > 0).sum(0)).ravel() >= 20].copy()
embed(iso_e, n_hvg=5000, n_pcs=30, resolution=0.5)

# isoform-usage embedding
u = usage_features(iso, min_gene_cells=0.2, min_iso_frac=0.05, max_genes=3000)
embed_dense(u, n_pcs=30, resolution=0.5)
print("usage features", u.shape, "genes", u.var.gene.nunique())

emb = []
for name, a in [("gene expression", g), ("isoform expression", iso_e), ("isoform usage", u)]:
    d = pd.DataFrame(a.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=a.obs_names)
    d["leiden"] = a.obs.leiden.values; d["sample"] = a.obs["sample"].values; d["space"] = name
    emb.append(d)
emb = pd.concat(emb)
emb = emb.join(g.obs[["phase", "S_score", "G2M_score", "total_counts", "pct_counts_mt"]])
emb.to_csv(f"{OUT}/tables/umap_coords.tsv", sep="\t")
iso_e.obs[["leiden"]].rename(columns={"leiden": "leiden_iso"}).join(u.obs[["leiden"]].rename(columns={"leiden": "leiden_usage"})) \
    .join(g.obs[["leiden"]]).to_csv(f"{OUT}/tables/clusters.tsv", sep="\t")
emb = emb.sample(frac=1, random_state=0)
for col in ["sample", "leiden", "phase"]:
    p = (ggplot(emb, aes("UMAP1", "UMAP2", color=col)) + geom_point(size=0.1, alpha=0.5) + facet_wrap("~space", scales="free")
         + theme_bw() + guides(color=guide_legend(override_aes={"size": 3, "alpha": 1})))
    p.save(f"{FIG}/umap_{col}.png", width=13, height=4, dpi=150)

# osteoblast / proliferation / temperature marker genes on the gene UMAP
markers = ["RUNX2", "SP7", "ALPL", "COL1A1", "COL1A2", "SPP1", "BGLAP", "IBSP", "SPARC", "POSTN", "MKI67", "TOP2A",
           "HSPA1A", "DNAJB1", "RBM3", "CIRBP", "CDKN1A", "MDM2"]
markers = [x for x in markers if x in g.var_names]
d = pd.DataFrame(g.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=g.obs_names)
ex = sc.get.obs_df(g, keys=markers)
d = d.join(ex).melt(id_vars=["UMAP1", "UMAP2"], var_name="gene", value_name="log_expr")
p = (ggplot(d.sort_values("log_expr"), aes("UMAP1", "UMAP2", color="log_expr")) + geom_point(size=0.05)
     + facet_wrap("~gene", ncol=6) + scale_color_cmap("viridis") + theme_bw())
p.save(f"{FIG}/umap_markers.png", width=15, height=8, dpi=130)

# cluster composition and contingency between spaces
ct = pd.crosstab(g.obs.leiden, g.obs["sample"]); ct.to_csv(f"{OUT}/tables/leiden_by_sample.tsv", sep="\t"); print(ct)
print(pd.crosstab(g.obs["sample"], g.obs.phase, normalize="index"))
