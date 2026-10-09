"""Alternative events (SUPPA2 local events, TSS and polyA-site clusters), detectability and day4 vs day0 tests."""
import sys, os, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from plotnine import *
from hfob_kinnex.annot import read_gtf_exons
from hfob_kinnex.events import (transcript_structure, end_clusters, classify_terminal, read_ioe, incidence,
                                group_incidence)
from hfob_kinnex.stats import quasibinomial_two_group

J, OUT, SUPPA = sys.argv[1], sys.argv[2], sys.argv[3]
EV = f"{J}/events"; os.makedirs(EV, exist_ok=True)
FIG, TAB = f"{OUT}/figures", f"{OUT}/tables"
MIN_READS_ISO = 10      # pooled reads for an isoform to enter event definition
DPSI = 0.1

iso = ad.read_h5ad(f"{J}/isoforms.h5ad")
qc = pd.read_csv(f"{TAB}/cell_qc.tsv", sep="\t", index_col=0)
iso = iso[qc.index[qc.qc_pass]].copy()
X = iso.X.tocsr()
day4 = (iso.obs["sample"] == "day4").values
tot = np.asarray(X.sum(0)).ravel()
expressed = iso.var_names[tot >= MIN_READS_ISO]
gene = iso.var.gene.astype(str)
print("cells", iso.n_obs, "isoforms", iso.n_vars, "expressed isoforms", len(expressed), flush=True)

# ---- SUPPA2 local events on expressed isoforms
gtf = f"{EV}/expressed.gtf"
if not os.path.exists(gtf):
    keep = set(expressed)
    with open(f"{J}/joint.sorted.filtered_lite.gff") as f, open(gtf, "w") as o:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            tid = p[8].split('transcript_id "')[1].split('"')[0]
            if tid in keep and p[2] == "exon":
                o.write(line)
if not os.path.exists(f"{EV}/suppa_SE_strict.ioe"):
    subprocess.run([sys.executable, f"{SUPPA}/suppa.py", "generateEvents", "-i", gtf, "-o", f"{EV}/suppa",
                    "-f", "ioe", "-e", "SE", "SS", "MX", "RI", "FL"], check=True)
ioe = pd.concat([read_ioe(f"{EV}/suppa_{t}_strict.ioe") for t in ["SE", "A5", "A3", "MX", "RI", "AF", "AL"]],
                ignore_index=True)
# map PB gene id -> associated gene name via the inclusion transcripts
ioe["gene_name"] = [gene.get(l[0], "NA") for l in ioe.total]
Inc, Tot = incidence(ioe.incl, iso.var_names), incidence(ioe.total, iso.var_names)
# restrict to events with pooled support (>= 50 reads in each day) before forming cell x event matrices
# per-day pooled inclusion / total reads for all events (used for detection counts)
ioe_all = ioe.copy()
py = np.vstack([np.asarray(X[m].sum(0)).ravel() @ Inc for m in (~day4, day4)])
pn = np.vstack([np.asarray(X[m].sum(0)).ravel() @ Tot for m in (~day4, day4)])
keep_ev = (pn >= 50).all(0)
print("SUPPA events", len(ioe), "with >=50 reads per day", keep_ev.sum(), flush=True)
ioe = ioe.loc[keep_ev].reset_index(drop=True)
Inc, Tot = Inc[:, keep_ev], Tot[:, keep_ev]
Y, N = (X @ Inc).tocsc(), (X @ Tot).tocsc()

# ---- TSS / polyA-site clusters on expressed isoforms
ex = read_gtf_exons(gtf)
tx = transcript_structure(ex)
cnt = pd.Series(tot, index=iso.var_names)
res_ends = {}
for which, win in [("tss", 50), ("tes", 100)]:
    lab, cl = end_clusters(tx, gene, cnt, which=which, window=win)
    # end support: TSS within a refTSS (CAGE) peak; 3' end with an upstream polyA signal
    sup = (iso.var.within_CAGE_peak.astype(str) == "True") if which == "tss" else iso.var.polyA_motif.astype(str).ne("nan")
    sw = pd.DataFrame({"cluster": lab, "w": cnt.loc[lab.index] * sup.loc[lab.index].astype(float)}).groupby("cluster").w.sum()
    cl["supported_frac"] = (sw.reindex(cl.index) / cl["count"].replace(0, np.nan)).values
    term = classify_terminal(cl, min_frac=0.1, min_count=20)
    cl.to_csv(f"{TAB}/{which}_clusters.tsv", sep="\t")
    term.to_csv(f"{TAB}/{which}_alt_genes.tsv", sep="\t", index=False)
    M, groups = group_incidence(lab, iso.var_names)
    G, gnames = group_incidence(gene.loc[lab.index], iso.var_names)
    gidx = gnames.get_indexer(cl.loc[groups, "gene"])
    Yc = (X @ M).tocsr()
    Nc = ((X @ G).tocsc()[:, gidx]).tocsr()
    res_ends[which] = (Yc, Nc, cl.loc[groups].rename_axis("cluster").reset_index(), term)

# ---- detectability: pooled and per day. An event is detected if total reads >= 20 and 0.05 <= PSI <= 0.95
def detected(Y, N, mask, min_n=20, lo=0.05, hi=0.95):
    y = np.asarray(Y[mask].sum(0)).ravel(); n = np.asarray(N[mask].sum(0)).ravel()
    psi = np.divide(y, n, out=np.full_like(y, np.nan, dtype=float), where=n > 0)
    return (n >= min_n) & (psi >= lo) & (psi <= hi)

rows = []
allc = np.ones(iso.n_obs, bool)
def detected_sums(y, n, min_n=20, lo=0.05, hi=0.95):
    psi = np.divide(y, n, out=np.full_like(y, np.nan, dtype=float), where=n > 0)
    return (n >= min_n) & (psi >= lo) & (psi <= hi)

for t in ["SE", "A5", "A3", "MX", "RI", "AF", "AL"]:
    m = (ioe_all.type == t).values
    for s, (y, n) in [("pooled", (py.sum(0), pn.sum(0))), ("day0", (py[0], pn[0])), ("day4", (py[1], pn[1]))]:
        d = detected_sums(y[m], n[m])
        rows.append((t, s, m.sum(), d.sum(), ioe_all.loc[m].loc[d].gene_name.nunique()))
for which, label in [("tss", "alt TSS"), ("tes", "alt polyA")]:
    Yc, Nc, cl, term = res_ends[which]
    for s, mask in [("pooled", allc), ("day0", ~day4), ("day4", day4)]:
        d = detected(Yc, Nc, mask, lo=0.1, hi=0.9)
        # genes with >= 2 detected clusters
        gd = cl.loc[d].groupby("gene").size()
        rows.append((label, s, cl.gene.nunique(), (gd >= 2).sum(), (gd >= 2).sum()))
        if s == "pooled":
            gs = cl.loc[d & (cl.supported_frac >= 0.5).values].groupby("gene").size()
            sl = "CAGE-supported" if which == "tss" else "polyA-signal-supported"
            rows.append((f"{label} ({sl})", s, cl.gene.nunique(), (gs >= 2).sum(), (gs >= 2).sum()))
    rows.append((label + " (distinct terminal exon)", "pooled", cl.gene.nunique(), term.distinct_exons.sum(), term.distinct_exons.sum()))
    rows.append((label + " (tandem, same exon)", "pooled", cl.gene.nunique(), term.tandem.sum(), term.tandem.sum()))
det = pd.DataFrame(rows, columns=["event_type", "cells", "n_defined", "n_detected", "n_genes_detected"])
det.to_csv(f"{TAB}/event_detection.tsv", sep="\t", index=False)
print(det.to_string(), flush=True)

# ---- tests day4 vs day0
r = quasibinomial_two_group(Y, N, day4)
r = pd.concat([ioe[["event_id", "type", "gene_name"]], r], axis=1)
r["sig"] = (r.padj < 0.05) & (r.dPSI.abs() >= DPSI)
r.to_csv(f"{TAB}/diff_events.tsv.gz", sep="\t", index=False)
summ = [("local:" + t, d.testable.sum(), d.sig.sum(), d[d.sig].gene_name.nunique()) for t, d in r.groupby("type")]
for which, label in [("tss", "TSS cluster"), ("tes", "polyA cluster")]:
    Yc, Nc, cl, term = res_ends[which]
    multi = cl.groupby("gene")["gene"].transform("size").values >= 2
    rr = quasibinomial_two_group(Yc[:, multi], Nc[:, multi], day4)
    rr = pd.concat([cl.loc[multi].reset_index(drop=True), rr], axis=1)
    rr["sig"] = (rr.padj < 0.05) & (rr.dPSI.abs() >= DPSI)
    rr.to_csv(f"{TAB}/diff_{which}.tsv.gz", sep="\t", index=False)
    summ.append((label, rr.testable.sum(), rr.sig.sum(), rr[rr.sig].gene.nunique()))

# isoform-level DTU with the same model (isoform count vs gene total)
G, gnames = group_incidence(gene, iso.var_names)
niso = gene.map(gene[tot >= MIN_READS_ISO].value_counts()).fillna(0).values
multi = (niso >= 2) & (tot >= MIN_READS_ISO)
gidx = gnames.get_indexer(gene[multi])
Ng = (X @ G).tocsc()[:, gidx]
ri = quasibinomial_two_group(X.tocsc()[:, multi], Ng, day4)
ri.insert(0, "isoform", iso.var_names[multi]); ri.insert(1, "gene", gene.values[multi])
ri = pd.concat([ri, iso.var.loc[iso.var_names[multi], ["structural_category", "associated_transcript", "subcategory",
                                                     "coding", "predicted_NMD"]].reset_index(drop=True)], axis=1)
ri["sig"] = (ri.padj < 0.05) & (ri.dPSI.abs() >= DPSI)
ri.to_csv(f"{TAB}/diff_isoform_usage_qb.tsv.gz", sep="\t", index=False)
summ.append(("isoform usage (DTU)", ri.testable.sum(), ri.sig.sum(), ri[ri.sig].gene.nunique()))
summ = pd.DataFrame(summ, columns=["test", "n_tested", "n_sig", "n_genes_sig"])
summ.to_csv(f"{TAB}/diff_splicing_summary.tsv", sep="\t", index=False)
print(summ.to_string())

# volcano-style plot of dPSI vs significance for local events
d = r[r.testable].copy(); d["mlog10p"] = -np.log10(d.padj.clip(lower=1e-300))
p = (ggplot(d, aes("dPSI", "mlog10p", color="sig")) + geom_point(size=0.4, alpha=0.5) + facet_wrap("~type", nrow=2)
     + scale_color_manual(values={True: "#c0392b", False: "#999999"}) + theme_bw()
     + labs(x="ΔPSI (day4 − day0)", y="−log10 FDR"))
p.save(f"{FIG}/volcano_events.png", width=12, height=5.5, dpi=150)
