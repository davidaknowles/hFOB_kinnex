"""Gene-level differential expression day4 vs day0 (cell-level Wilcoxon + pseudobulk fold change) and pathway GSEA."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, gseapy as gp
from plotnine import *

J, OUT = sys.argv[1], sys.argv[2]
FIG, TAB = f"{OUT}/figures", f"{OUT}/tables"
g = ad.read_h5ad(f"{J}/genes.processed.h5ad")

# pseudobulk log2 CPM per day
C = g.layers["counts"]
pb = pd.DataFrame({s: np.asarray(C[(g.obs["sample"] == s).values].sum(0)).ravel() for s in ["day0", "day4"]},
                  index=g.var_names)
cpm = pb / pb.sum() * 1e6
pb_lfc = np.log2(cpm.day4 + 1) - np.log2(cpm.day0 + 1)
pct = pd.DataFrame({s: np.asarray((C[(g.obs["sample"] == s).values] > 0).mean(0)).ravel() for s in ["day0", "day4"]},
                   index=g.var_names)

sc.tl.rank_genes_groups(g, "sample", groups=["day4"], reference="day0", method="wilcoxon", tie_correct=True)
de = sc.get.rank_genes_groups_df(g, group="day4").set_index("names")
de = de.join(pd.DataFrame({"cpm_day0": cpm.day0, "cpm_day4": cpm.day4, "pb_log2FC": pb_lfc,
                           "pct_day0": pct.day0, "pct_day4": pct.day4}))
de["sig"] = (de.pvals_adj < 0.05) & (de.pb_log2FC.abs() >= 0.5) & (de[["pct_day0", "pct_day4"]].max(1) >= 0.1)
de["direction"] = np.where(~de.sig, "ns", np.where(de.pb_log2FC > 0, "up day4", "down day4"))
de.sort_values("scores", ascending=False).to_csv(f"{TAB}/de_genes.tsv.gz", sep="\t")
print(de.direction.value_counts())

lab = ["RUNX2", "SP7", "ALPL", "COL1A1", "SPP1", "IBSP", "BGLAP", "MKI67", "TOP2A", "HSPA1A", "HSPA1B", "DNAJB1",
       "RBM3", "CIRBP", "CDKN1A", "MDM2", "POSTN", "DKK1", "MGP", "TNFRSF11B"]
d = de.copy(); d["mlog10p"] = -np.log10(d.pvals_adj.clip(lower=1e-300)); d["gene"] = d.index
hi = d[d.sig & (d.cpm_day0 + d.cpm_day4 > 100)]
top = set(hi.nlargest(8, "pb_log2FC").index) | set(hi.nsmallest(8, "pb_log2FC").index)
lab = ["RUNX2", "ALPL", "SPP1", "BGLAP", "COL1A2", "POSTN", "MKI67", "CCNB1", "RBM3", "CDKN1A", "MDM2"]
d["label"] = np.where(d.index.isin(lab) | d.index.isin(top), d.index, "")
p = (ggplot(d, aes("pb_log2FC", "mlog10p", color="direction")) + geom_point(size=0.4, alpha=0.5)
     + geom_text(d[d.label != ""], aes(label="label"), size=7, color="black",
                 adjust_text={"arrowprops": {"arrowstyle": "-", "color": "grey", "lw": 0.5}})
     + scale_color_manual(values={"up day4": "#c0392b", "down day4": "#2e86c1", "ns": "#aaaaaa"})
     + theme_bw() + labs(x="pseudobulk log2FC (day4 / day0)", y="−log10 FDR (Wilcoxon)"))
p.save(f"{FIG}/volcano_genes.png", width=8, height=6, dpi=150)

# marker panel by day
panels = {"osteoblast": ["RUNX2", "SP7", "ALPL", "COL1A1", "COL1A2", "SPP1", "IBSP", "BGLAP", "SPARC", "POSTN", "TNFRSF11B", "MGP"],
          "proliferation": ["MKI67", "TOP2A", "CCNB1", "PCNA", "MCM2", "E2F1"],
          "heat shock": ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "HSP90AA1", "HSPB1"],
          "cold shock": ["RBM3", "CIRBP"],
          "p53 targets": ["CDKN1A", "MDM2", "GDF15", "BAX", "TP53I3", "SESN1"]}
rows = [(k, x) for k, v in panels.items() for x in v if x in de.index]
pm = pd.DataFrame(rows, columns=["panel", "gene"]).join(de[["pb_log2FC", "pvals_adj", "cpm_day0", "cpm_day4", "pct_day0", "pct_day4"]], on="gene")
pm.to_csv(f"{TAB}/marker_panel.tsv", sep="\t", index=False)
print(pm.to_string())
pm["gene"] = pd.Categorical(pm.gene, categories=pm.gene[::-1])
p = (ggplot(pm, aes("pb_log2FC", "gene", fill="panel")) + geom_col() + geom_vline(xintercept=0)
     + theme_bw() + labs(x="pseudobulk log2FC (day4 / day0)", y=""))
p.save(f"{FIG}/marker_panel.png", width=6, height=8, dpi=150)

# GSEA on pseudobulk fold change of genes expressed in >= 5% of cells in either day
rnk = de[(de[["pct_day0", "pct_day4"]].max(1) >= 0.05) & ~de.index.str.startswith(("ENSG", "novel"))].pb_log2FC
rnk = rnk.sort_values(ascending=False)
out = []
for lib in ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023"]:
    gs = gp.get_library(lib, organism="Human")
    r = gp.prerank(rnk=rnk, gene_sets=gs, min_size=15, max_size=500, permutation_num=2000, seed=0, threads=8,
                   outdir=None, verbose=False).res2d
    r["library"] = lib
    out.append(r)
gsea = pd.concat(out)
gsea["FDR q-val"] = gsea["FDR q-val"].astype(float); gsea["NES"] = gsea["NES"].astype(float)
gsea.sort_values("NES").to_csv(f"{TAB}/gsea_prerank.tsv", sep="\t", index=False)
top = gsea[(gsea.library == "MSigDB_Hallmark_2020") & (gsea["FDR q-val"] < 0.05)].sort_values("NES")
print(top[["Term", "NES", "FDR q-val"]].to_string())
top["Term"] = pd.Categorical(top.Term, categories=top.Term)
p = (ggplot(top, aes("NES", "Term", fill="NES > 0")) + geom_col() + theme_bw() + labs(y="", x="NES (day4 vs day0)")
     + scale_fill_manual(values={True: "#c0392b", False: "#2e86c1"}) + theme(legend_position="none"))
p.save(f"{FIG}/gsea_hallmark.png", width=7, height=max(3, 0.25 * len(top)), dpi=150)
