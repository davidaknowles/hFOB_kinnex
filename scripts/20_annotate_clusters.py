"""Annotate gene-expression Leiden clusters: composition, QC, cell-cycle phase, curated state scores, markers, Hallmark ORA."""
import sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, gseapy as gp
from plotnine import *

J, OUT = sys.argv[1], sys.argv[2]
FIG, TAB = f"{OUT}/figures", f"{OUT}/tables"
g = ad.read_h5ad(f"{J}/genes.processed.h5ad")
g.obs["leiden"] = g.obs.leiden.astype(str)
ribo = g.var_names.str.match(r"^RP[LS]\d") & ~g.var_names.str.contains("_")
g.obs["pct_ribo"] = np.asarray(g.layers["counts"][:, ribo].sum(1)).ravel() / g.obs.total_counts * 100

STATES = {
    "S phase": None, "G2M": None,
    "interferon": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "IFI6", "ISG20", "OASL", "IFI44", "RSAD2", "STAT1"],
    "immediate early / stress": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "IER2", "DUSP1", "KLF6"],
    "DNA damage / p53": ["CDKN1A", "MDM2", "GDF15", "TP53I3", "SESN1", "BAX", "GADD45A", "FDXR", "RRM2B", "PLK3"],
    "SASP": ["IL6", "CXCL8", "CXCL1", "CXCL2", "SERPINE1", "IL1B", "MMP1", "MMP3", "IGFBP3", "SAA1"],
    "myofibroblast": ["ACTA2", "TAGLN", "MYL9", "CNN1", "TPM2", "CALD1"],
    "ECM": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "POSTN", "SPARC", "COL4A1", "COL4A2"],
    "osteoblast": ["RUNX2", "SP7", "ALPL", "SPP1", "IBSP", "BGLAP", "COL1A1", "TNFRSF11B", "DLX5"],
    "hypoxia / glycolysis": ["VEGFA", "BNIP3", "CA9", "PGK1", "LDHA", "ENO1", "SLC2A1", "NDRG1", "P4HA1"],
    "histone (replication)": ["H4C3", "H1-2", "H1-4", "H2AC14", "H1-5", "H2BC12", "H3C2"],
    "UPR / ER stress": ["HSPA5", "DDIT3", "XBP1", "ATF4", "HERPUD1", "SEL1L", "DNAJB9"],
    "heat shock": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "HSPA6", "BAG3"],
}
for k, v in STATES.items():
    if v is None:
        continue
    v = [x for x in v if x in g.var_names]
    sc.tl.score_genes(g, v, score_name=k, random_state=0)
score_cols = ["S_score", "G2M_score"] + [k for k, v in STATES.items() if v is not None]

grp = g.obs.groupby("leiden")
summ = grp.agg(n=("sample", "size"), frac_day4=("sample", lambda x: (x == "day4").mean()),
               umis=("total_counts", "median"), genes=("n_genes_by_counts", "median"),
               pct_mt=("pct_counts_mt", "median"), pct_ribo=("pct_ribo", "median"))
ph = pd.crosstab(g.obs.leiden, g.obs.phase, normalize="index").add_prefix("frac_")
summ = summ.join(ph).join(grp[score_cols].mean())
summ = summ[summ.n >= 50]

# markers: each cluster vs the other clusters of its majority day (cluster 2 is mostly day4)
rows, ora = [], []
hall = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
bg = list(g.var_names[np.asarray((g.layers["counts"] > 0).mean(0)).ravel() >= 0.05])
for cl in summ.index:
    day = "day4" if summ.loc[cl, "frac_day4"] > 0.5 else "day0"
    a = g[(g.obs["sample"] == day).values].copy()
    a.obs["grp"] = np.where(a.obs.leiden == cl, "in", "out")
    sc.tl.rank_genes_groups(a, "grp", groups=["in"], reference="out", method="wilcoxon", pts=True)
    d = sc.get.rank_genes_groups_df(a, group="in")
    d = d[~d.names.str.match(r"^(RP[LS]\d|MT-|novelGene|ENSG)") & ~d.names.str.contains("_")]
    up = d[(d.pvals_adj < 0.05) & (d.logfoldchanges > 0.5) & (d.pct_nz_group > 0.2)].sort_values("scores", ascending=False)
    rows.append((cl, day, ", ".join(up.names.head(20))))
    if len(up) >= 10:
        r = gp.enrich(gene_list=list(up.names.head(150)), gene_sets=hall, background=bg, outdir=None, verbose=False).results
        r = r[r["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").head(4)
        ora.append((cl, "; ".join(f"{t} ({p:.0e})" for t, p in zip(r.Term, r["Adjusted P-value"]))))
mk = pd.DataFrame(rows, columns=["leiden", "vs_day", "top_markers"]).set_index("leiden")
mk = mk.join(pd.DataFrame(ora, columns=["leiden", "hallmark_ORA"]).set_index("leiden"))
summ = summ.join(mk)
summ.to_csv(f"{TAB}/cluster_annotation.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
print(summ.drop(columns=["top_markers", "hallmark_ORA"]).round(2).to_string())
print(summ[["vs_day", "top_markers", "hallmark_ORA"]].to_string())

# heatmap of z-scored state scores per cluster
h = summ[score_cols]
h = ((h - h.mean()) / h.std()).reset_index().melt(id_vars="leiden", var_name="state", value_name="z")
order = summ.sort_values("frac_day4").index.tolist()
h["leiden"] = pd.Categorical(h.leiden, categories=order)
h["state"] = pd.Categorical(h.state, categories=score_cols[::-1])
lab = {c: f"{c} ({'d4' if summ.loc[c, 'frac_day4'] > 0.5 else 'd0'}, n={summ.loc[c, 'n']})" for c in order}
p = (ggplot(h, aes("leiden", "state", fill="z")) + geom_tile() + scale_fill_gradient2(low="#2e86c1", mid="white", high="#c0392b")
     + scale_x_discrete(labels=[lab[c] for c in order]) + theme_bw() + theme(axis_text_x=element_text(rotation=45, ha="right"))
     + labs(x="Leiden cluster", y="", fill="z (across\nclusters)"))
p.save(f"{FIG}/cluster_state_scores.png", width=8, height=5.5, dpi=150)
d = pd.DataFrame(g.obsm["X_umap"], columns=["UMAP1", "UMAP2"]).assign(leiden=g.obs.leiden.values)
cent = d.groupby("leiden")[["UMAP1", "UMAP2"]].median().reset_index()
p = (ggplot(d.sample(frac=1, random_state=0), aes("UMAP1", "UMAP2", color="leiden")) + geom_point(size=0.1, alpha=0.5)
     + geom_text(cent, aes(label="leiden"), color="black", size=11) + theme_bw()
     + guides(color=guide_legend(override_aes={"size": 3, "alpha": 1})))
p.save(f"{FIG}/umap_leiden_labeled.png", width=6.5, height=5, dpi=150)
