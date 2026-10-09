"""Per-cell and kNN-smoothed PSI on the gene-expression UMAP for example events with known biology."""
import sys, os, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, anndata as ad
from plotnine import *
from hfob_kinnex.annot import read_gtf_exons
from hfob_kinnex.events import read_ioe, transcript_structure, end_clusters, knn_smooth_ratio

J, OUT = sys.argv[1], sys.argv[2]
EX = [  # label, kind, id ; kind is a SUPPA event id or a 3' end cluster (gene|tesK)
    ("RPS24 alt exon inclusion", "suppa", "PB.190333;SE:chr10:78037304-78040204:78040225-78040615:+"),
    ("CIRBP distal last exon", "tes", "CIRBP|tes9"),
    ("CIRBP 3' intron retention", "suppa", "PB.275388;RI:chr19:1274307:1274440-1274649:1274807:+"),
    ("H1-2 polyadenylated 3' end", "tes", "H1-2|tes1"),
    ("H2AC6 proximal 3' end", "tes", "H2AC6|tes0"),
    ("TPM2 exon 6a vs 6b", "suppa", "PB.172564;MX:chr9:35684550-35684732:35684807-35685269:35684550-35685064:35685139-35685269:-"),
    ("TPM1 alt first exon", "suppa", "PB.249466;AF:chr15:63042481:63042943-63056985:63048209:63048707-63056985:+"),
    ("CALD1 alt first exon", "suppa", "PB.151068;AF:chr7:134646856:134646937-134843884:134779603:134779749-134843884:+"),
    ("HNRNPDL exon inclusion", "suppa", "PB.81987;SE:chr4:82424883-82425563:82425667-82426037:-"),
    ("SRSF5 alt 5' splice site", "suppa", "PB.241751;A5:chr14:69767380-69768138:69767251-69768138:+"),
    ("SRSF2 poison exon", "suppa", "PB.266522;SE:chr17:76735158-76735772:76735875-76736154:-"),
    ("SRSF7 intron retention", "suppa", "PB.28044;RI:chr2:38746694:38746747-38748056:38748157:-"),
]

g = ad.read_h5ad(f"{J}/genes.processed.h5ad")
iso = ad.read_h5ad(f"{J}/isoforms.h5ad")[g.obs_names]
X = iso.X.tocsc()
gene = iso.var.gene.astype(str)
ioe = pd.concat([read_ioe(f) for f in glob.glob(f"{J}/events/suppa_*_strict.ioe")]).set_index("event_id")
ex = read_gtf_exons(f"{J}/events/expressed.gtf")
cnt = pd.Series(np.asarray(X.sum(0)).ravel(), index=iso.var_names)

def cols(ids):
    j = iso.var_names.get_indexer(ids)
    return np.asarray(X[:, j[j >= 0]].sum(1)).ravel()

rows, pb = [], []
for label, kind, eid in EX:
    if kind == "suppa":
        e = ioe.loc[eid]
        y, n = cols(e.incl), cols(e.total)
    else:
        gname = eid.split("|")[0]
        ids = gene.index[gene == gname]
        tx = transcript_structure(ex[ex.transcript_id.isin(ids)])
        lab, _ = end_clusters(tx, gene, cnt, which="tes", window=100)
        y, n = cols(lab.index[lab == eid]), cols(lab.index)
    psi_s, ns = knn_smooth_ratio(g.obsp["connectivities"], y, n)
    rows.append(pd.DataFrame({"event": label, "UMAP1": g.obsm["X_umap"][:, 0], "UMAP2": g.obsm["X_umap"][:, 1],
                              "sample": g.obs["sample"].values, "leiden": g.obs.leiden.values,
                              "y": y, "n": n, "psi_smooth": psi_s}, index=g.obs_names))
    for s in ["day0", "day4"]:
        m = (g.obs["sample"] == s).values
        pb.append((label, eid, s, y[m].sum(), n[m].sum(), y[m].sum() / max(n[m].sum(), 1), (n[m] > 0).mean()))
d = pd.concat(rows)
d["event"] = pd.Categorical(d.event, categories=[e[0] for e in EX])
pb = pd.DataFrame(pb, columns=["event", "event_id", "sample", "inclusion_reads", "total_reads", "PSI", "frac_cells_covered"])
pb.to_csv(f"{OUT}/tables/known_biology_psi.tsv", sep="\t", index=False)
print(pb.round(3).to_string())

FIG = f"{OUT}/figures"
dd = d.dropna(subset=["psi_smooth"]).sample(frac=1, random_state=0)
p = (ggplot(dd, aes("UMAP1", "UMAP2", color="psi_smooth")) + geom_point(size=0.08) + facet_wrap("~event", ncol=4)
     + scale_color_cmap("viridis", limits=(0, 1), name="PSI\n(kNN-smoothed)") + theme_bw()
     + theme(strip_text=element_text(size=8)))
p.save(f"{FIG}/umap_psi_known_biology.png", width=14, height=10, dpi=130)
# per-cell PSI distribution (cells with >= 3 event reads) with pseudobulk PSI
c = d[d.n >= 3].assign(psi=lambda x: x.y / x.n)
p = (ggplot(c, aes("sample", "psi", fill="sample")) + geom_violin(scale="width") + facet_wrap("~event", ncol=6)
     + geom_point(pb.assign(event=pd.Categorical(pb.event, categories=[e[0] for e in EX])), aes("sample", "PSI"),
                  inherit_aes=False, color="black", size=2.5, shape="D")
     + theme_bw() + theme(legend_position="none", strip_text=element_text(size=7)) + labs(x="", y="per-cell PSI (≥3 reads); ◆ pseudobulk"))
p.save(f"{FIG}/psi_violin_known_biology.png", width=14, height=5.5, dpi=130)
# per-cluster pseudobulk PSI
cl = d.groupby(["event", "sample", "leiden"], observed=True)[["y", "n"]].sum().reset_index()
cl = cl[cl.n >= 20].assign(PSI=lambda x: x.y / x.n)
p = (ggplot(cl, aes("leiden", "PSI", color="sample", size="n")) + geom_point() + facet_wrap("~event", ncol=4, scales="free_y")
     + scale_size_continuous(range=(0.8, 4), name="reads") + theme_bw() + theme(strip_text=element_text(size=8))
     + labs(x="Leiden cluster (gene expression)", y="pseudobulk PSI"))
p.save(f"{FIG}/psi_by_cluster_known_biology.png", width=14, height=9, dpi=130)
