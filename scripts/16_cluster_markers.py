"""Marker genes per Leiden cluster (Wilcoxon, one-vs-rest within each day) and osteoblast/p53/cycle scores."""
import sys
import pandas as pd, scanpy as sc, anndata as ad
J, OUT = sys.argv[1], sys.argv[2]
g = ad.read_h5ad(f"{J}/genes.processed.h5ad")
rows = []
for s in ["day0", "day4"]:
    a = g[g.obs["sample"] == s].copy()
    a = a[:, a.var_names].copy()
    keep = a.obs.leiden.value_counts(); keep = keep[keep >= 50].index
    a = a[a.obs.leiden.isin(keep)].copy(); a.obs.leiden = a.obs.leiden.astype(str)
    sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon")
    for cl in keep:
        d = sc.get.rank_genes_groups_df(a, group=str(cl)).head(15)
        rows.append((s, cl, (a.obs.leiden == str(cl)).sum(), ", ".join(d.names)))
m = pd.DataFrame(rows, columns=["sample", "leiden", "n_cells", "top_markers"])
m.to_csv(f"{OUT}/tables/cluster_markers.tsv", sep="\t", index=False)
pd.set_option("display.max_colwidth", 200); pd.set_option("display.width", 250)
print(m.to_string())
sets = {"osteoblast": ["ALPL", "SPP1", "COL1A1", "COL1A2", "BGLAP", "IBSP", "RUNX2", "SP7", "POSTN"],
        "p53": ["CDKN1A", "MDM2", "GDF15", "TP53I3", "SESN1", "BAX"]}
for k, v in sets.items():
    sc.tl.score_genes(g, [x for x in v if x in g.var_names], score_name=k)
print(g.obs.groupby(["sample", "leiden"], observed=True)[["osteoblast", "p53", "S_score", "G2M_score"]].mean().round(2))
