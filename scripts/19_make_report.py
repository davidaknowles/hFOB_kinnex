"""Self-contained HTML report (figures embedded as base64) from results/figures and results/tables."""
import sys, os, base64, html
import numpy as np
import pandas as pd

R = sys.argv[1]
OUT = f"{R}/report.html"
FIG, TAB = f"{R}/figures", f"{R}/tables"


MODE = "html"  # "html": self-contained, images embedded; "md": GitHub-rendered markdown, images linked


def img(path, caption, width="100%"):
    if MODE == "md":
        return f'<p><img src="figures/{path}" width="{width}" alt="{html.escape(caption)}"><br><em>{caption}</em></p>\n'
    with open(f"{FIG}/{path}", "rb") as f:
        b = base64.b64encode(f.read()).decode()
    return (f'<figure><img src="data:image/png;base64,{b}" style="max-width:{width}" alt="{html.escape(caption)}">'
            f"<figcaption>{caption}</figcaption></figure>")


def table(df, digits=3, index=False):
    return df.round(digits).to_html(index=index, border=0, classes="tbl", na_rep="")


qc = pd.read_csv(f"{TAB}/cell_qc_summary.tsv", sep="\t")
rl = pd.read_csv(f"{TAB}/read_length_summary.tsv", sep="\t")[["type", "sample", "10%", "50%", "90%"]]
det = pd.read_csv(f"{TAB}/event_detection.tsv", sep="\t")
det = det[det.cells == "pooled"].drop(columns="cells")
ds = pd.read_csv(f"{TAB}/diff_splicing_summary.tsv", sep="\t")
mk = pd.read_csv(f"{TAB}/marker_panel.tsv", sep="\t")
gs = pd.read_csv(f"{TAB}/gsea_prerank.tsv", sep="\t")
gs = gs[(gs.library == "MSigDB_Hallmark_2020") & (gs["FDR q-val"] < 0.05)].sort_values("NES", ascending=False)[["Term", "NES", "FDR q-val"]]
kb = pd.read_csv(f"{TAB}/known_biology_psi.tsv", sep="\t")
kb = kb.pivot_table(index="event", columns="sample", values=["PSI", "total_reads"], sort=False)
kb.columns = [f"{a} {b}" for a, b in kb.columns]
kb["ΔPSI"] = kb["PSI day4"] - kb["PSI day0"]
kb = kb[["PSI day0", "PSI day4", "ΔPSI", "total_reads day0", "total_reads day4"]].astype({"total_reads day0": int, "total_reads day4": int})
# cluster annotation
ca = pd.read_csv(f"{TAB}/cluster_annotation.tsv", sep="\t", dtype={"leiden": str}).set_index("leiden")
cl_lab = pd.read_csv(f"{TAB}/cluster_labels.tsv", sep="\t", dtype={"leiden": str}).set_index("leiden")
ca = cl_lab.join(ca)
ca["day"] = np.where(ca.frac_day4 > 0.5, "day 4", "day 0")
ca["markers (top 10)"] = ca.top_markers.str.split(", ").str[:10].str.join(", ")
cat = ca.reset_index()[["leiden", "day", "label", "n", "frac_day4", "umis", "genes", "frac_G1", "frac_S", "frac_G2M",
                        "markers (top 10)", "hallmark_ORA"]].sort_values(["day", "leiden"])
cat = cat.astype({"n": int, "umis": int, "genes": int})

# sc-SHC results (if available)
SC = sys.argv[2] if len(sys.argv) > 2 else None  # sc-SHC output directory
def scshc_html():
    if SC is None or not os.path.exists(f"{SC}/day0_testClusters.tsv"):
        return "<p class=\"muted\">sc-SHC results not yet available.</p>"
    out, rows = [], []
    for s in ["day0", "day4"]:
        if not os.path.exists(f"{SC}/{s}_testClusters.tsv"):
            continue
        t = pd.read_csv(f"{SC}/{s}_testClusters.tsv", sep="\t", dtype=str)
        merged = t.groupby("scshc_merged").leiden.apply(lambda x: "+".join(sorted(x.unique(), key=int)))
        row = [s, t.leiden.nunique(), t.scshc_merged.nunique(), "; ".join(merged.values)]
        if os.path.exists(f"{SC}/{s}_denovo.tsv"):
            d = pd.read_csv(f"{SC}/{s}_denovo.tsv", sep="\t", dtype=str)
            row.append(d.scshc_denovo.nunique())
            ct = pd.crosstab(d.leiden, d.scshc_denovo)
            ct = ct.loc[ct.sum(1) >= 20]
            out.append(f"<h3>{s}: de novo sc-SHC clusters vs Leiden</h3>" + ct.to_html(border=0, classes="tbl"))
        else:
            row.append("running")
        rows.append(row)
    tb = pd.DataFrame(rows, columns=["day", "Leiden clusters tested (≥20 cells)", "significant clusters after sc-SHC merging",
                                     "groups (merged Leiden clusters)", "de novo sc-SHC clusters"])
    return table(tb) + "".join(out)

iso_figs = ["RPS24", "CIRBP", "H1-2", "TPM2", "TPM1", "CALD1", "HNRNPDL", "SRSF2", "SRSF3"]

css = """
:root { --bg:#ffffff; --fg:#1d1d1f; --muted:#5f6368; --line:#e3e3e3; --accent:#2e6da4; --card:#f7f7f8; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#16181b; --fg:#e8e8e8; --muted:#a0a4a8; --line:#33373c; --accent:#7fb2e5; --card:#1f2226; } }
:root[data-theme="dark"] { --bg:#16181b; --fg:#e8e8e8; --muted:#a0a4a8; --line:#33373c; --accent:#7fb2e5; --card:#1f2226; }
body { background:var(--bg); color:var(--fg); font:15px/1.55 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin:0; }
main { max-width:1100px; margin:0 auto; padding:24px 16px 64px; }
h1 { font-size:26px; margin:0 0 4px; } h2 { font-size:20px; margin:40px 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px; }
h3 { font-size:16px; margin:24px 0 6px; } p, li { max-width:900px; } .muted { color:var(--muted); }
nav { background:var(--card); border:1px solid var(--line); border-radius:8px; padding:8px 14px; margin:16px 0; }
nav a { color:var(--accent); margin-right:14px; text-decoration:none; white-space:nowrap; }
figure { margin:12px 0 20px; } figure img { width:100%; height:auto; background:#fff; border:1px solid var(--line); border-radius:6px; }
figcaption { color:var(--muted); font-size:13px; margin-top:4px; max-width:900px; }
.grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(320px, 1fr)); gap:12px; }
.tbl { border-collapse:collapse; font-size:13px; margin:8px 0 16px; display:block; overflow-x:auto; }
.tbl th, .tbl td { padding:4px 10px; border-bottom:1px solid var(--line); text-align:right; }
.tbl th:first-child, .tbl td:first-child { text-align:left; }
.callout { background:var(--card); border-left:3px solid var(--accent); padding:10px 14px; border-radius:4px; max-width:900px; }
"""

def build():
    return f"""
<h1>hFOB 1.19 Kinnex single-cell long-read: day 4 vs day 0</h1>
<p class="muted">PacBio Kinnex (10x 3′ v3.x) on hFOB 1.19 osteoblasts, one library per time point. Generated 2026-10-09.</p>
<nav><a href="#summary">Summary</a><a href="#methods">Methods</a><a href="#qc">Read QC</a><a href="#cells">Cells</a>
<a href="#subpops">Subpopulations</a><a href="#de">Expression</a><a href="#events">Events</a><a href="#biology">Known biology</a><a href="#caveats">Caveats</a></nav>

<h2 id="summary">Summary</h2>
<ul>
<li>16,180 cells (10,258 day 0, 5,939 day 4) from SMRT Link; 15,718 pass QC. One joint isoform set across both days (isoseq collapse + pigeon): 1.28M filtered isoforms.</li>
<li>cDNA reads are short (median ~530 nt) and 40% of reads are removed by pigeon (mostly intra-priming). Coverage is slightly 5′-biased: molecules reach the TSS but many 3′ ends fall early.</li>
<li>Days separate completely in gene expression, isoform expression and isoform usage space.</li>
<li>Day 4 shows p53 reactivation (CDKN1A, MDM2, GDF15), proliferation arrest, a modest early osteoblast program (ALPL, SPP1, COL1A2, POSTN) and loss of the cold-inducible RBM3. This matches inactivation of the temperature-sensitive SV40 large T at 39.5 °C plus early differentiation.</li>
<li>Isoform-level changes with known biology: temperature-dependent CIRBP 3′ processing, SR-protein autoregulatory exons, histone mRNA polyadenylation at cell-cycle exit, tropomyosin / caldesmon isoform switches, and a near-complete RPS24 alternative-exon switch.</li>
</ul>

<h2 id="methods">Methods</h2>
<ol>
<li>Real-cell reads from each SMRT Link mapped BAM (pbmm2, hg38) were renamed with a sample prefix on read and cell barcode, merged, and collapsed jointly with <code>isoseq collapse</code>. Isoforms were classified with <code>pigeon classify</code> (GENCODE v39, refTSS v3.3 CAGE peaks, polyA motifs) and filtered with <code>pigeon filter</code> (intra-priming, RT switching, junction support).</li>
<li>Counts: each read is one UMI-deduplicated molecule; cell × isoform counts come from read → cell barcode and read → isoform; gene counts sum isoforms of each associated gene.</li>
<li>Cells: per-sample 3-MAD filters on log UMIs, log genes (low) and mitochondrial percent (high). scanpy normalisation, seurat_v3 HVGs, PCA, UMAP and Leiden in three spaces. Isoform usage features are per-cell isoform/gene proportions centred on the pooled usage, for multi-isoform genes detected in ≥ 20% of cells.</li>
<li>Events: SUPPA2 local events (SE, A5, A3, MX, RI, AF, AL) on isoforms with ≥ 10 reads; TSS (50 nt) and 3′ end (100 nt) clusters per gene. Detected = ≥ 20 pooled reads with PSI in [0.05, 0.95].</li>
<li>Differential usage: quasi-binomial GLM with day as the only covariate, cells as observations (closed form, Pearson dispersion floored at 1), BH FDR &lt; 0.05 and |ΔPSI| ≥ 0.1. satuRn was run on the isoform set as a check.</li>
<li>Gene DE: cell-level Wilcoxon with pseudobulk log2 CPM fold change (FDR &lt; 0.05, |log2FC| ≥ 0.5, detected in ≥ 10% of cells). GSEA prerank on pseudobulk log2FC (MSigDB Hallmark, GO BP).</li>
</ol>

<h2 id="qc">Read QC</h2>
<div class="grid">{img("knee.png", "Barcode-rank plots. Red: barcodes called as cells by SMRT Link.")}
{img("genebody_coverage.png", "RSeQC gene-body coverage on canonical protein-coding transcripts (2% of reads).")}</div>
{img("read_lengths.png", "Read lengths: HiFi array reads (Kinnex concatemers), segmented S-reads, and deduplicated mapped reads from real cells.")}
{table(rl)}
<div class="grid">{img("read_filter_reasons.png", "Per-read pigeon filter outcome (QC-pass cells).")}
{img("read_structural_category.png", "Structural category of passing reads.")}</div>
<p>Among passing reads assigned to known transcripts, 5′-fragments (missing 3′ exons) are 8.5% / 7.3% and 3′-fragments (missing 5′ exons) 3.0% (day 0 / day 4). The mild 5′-high coverage therefore reflects early 3′ ends (residual internal priming, proximal polyA sites upstream of long annotated UTRs) rather than RT truncation.</p>

<h2 id="cells">Cells and embeddings</h2>
{table(qc)}
{img("cell_qc_violin.png", "Per-cell UMIs, genes, isoforms and mitochondrial percent (counts after pigeon filtering).")}
{img("umap_sample.png", "UMAPs on gene expression, isoform expression and isoform usage, coloured by day.")}
{img("umap_phase.png", "Cell-cycle phase (Seurat 2019 S/G2M genes). Day 0: 62% G2M; day 4: 46% G1.")}
{img("umap_markers.png", "Marker genes on the gene-expression UMAP (log-normalised).")}

<h2 id="subpops">Subpopulations</h2>
<p>Leiden clusters on the gene-expression graph (all cells, resolution 0.5) were annotated with within-day one-vs-rest Wilcoxon markers (excluding ribosomal, mitochondrial and unannotated genes), Hallmark over-representation of the top 150 markers, cell-cycle phase, and curated state scores (scanpy score_genes).</p>
<div class="grid">{img("umap_leiden_labeled.png", "Leiden clusters on the gene-expression UMAP.")}
{img("cluster_state_scores.png", "Mean state scores per cluster, z-scored across clusters.")}</div>
{table(cat, digits=2)}
<p>Day 0 subclusters are mainly cell-cycle phases (G1/S, S/G2 with histone mRNAs, G2/M), plus a mesenchymal / EMT-high group, a small interferon-response group (213 cells) and a low-complexity ribosome/OXPHOS-high group. Day 4 has an ECM-producing majority (POSTN, FN1, COL4A1, INHBA), a mesenchymal group with the highest osteoblast score (SPARC, COL1A1), a low-complexity group matching day 0 cluster 4, and a small p53/stress group (GADD45A, PLK2, AREG) that includes 26 day 0 cells. Part of day 0 cluster 3 bridges toward day 4 cluster 7, suggesting a day 0 subset already resembling the day 4 mesenchymal state. The two low-complexity clusters (fewest genes and UMIs, short housekeeping transcripts) may be smaller or lower-quality cells rather than a distinct biological state.</p>
<h3>Are the clusters statistically supported? sc-SHC</h3>
<p>sc-SHC (Grabski, Street &amp; Irizarry 2023, Nat Methods) tests each split of a cluster hierarchy against a null of a single population (Gaussian-copula model of the counts fitted to the merged cells), controlling the family-wise error rate (α = 0.05). It was run per day on raw gene counts (2,500 features, 30 PCs): <code>testClusters</code> on the Leiden labels (clusters with ≥ 20 cells in that day), which merges clusters that are not significantly different, and <code>scSHC</code> for de novo significance-based clustering.</p>
{scshc_html()}

<h2 id="de">Differential expression, day 4 vs day 0</h2>
<p>1,501 genes up and 1,553 down.</p>
<div class="grid">{img("volcano_genes.png", "Volcano: pseudobulk log2FC vs cell-level Wilcoxon FDR (capped at 1e-300).")}
{img("marker_panel.png", "Marker panels: osteoblast, proliferation, heat shock, cold shock, p53 targets.")}</div>
{table(mk)}
<h3>Hallmark GSEA (FDR &lt; 0.05)</h3>
{table(gs)}

<h2 id="events">Alternative isoform events</h2>
<h3>Detected events (pooled cells)</h3>
{table(det)}
<p>Requiring every TSS cluster to lie in a refTSS peak reduces alternative-TSS genes from 8,185 to 895, so most raw alternative TSS calls are likely 5′-truncated molecules. SUPPA AF/AL counts are inflated by combinatorial isoform pairs; the cluster-based counts are preferred.</p>
<h3>Differential usage (FDR &lt; 0.05, |ΔPSI| ≥ 0.1)</h3>
{table(ds)}
{img("volcano_events.png", "ΔPSI vs FDR for SUPPA local events.")}
<div class="grid">{img("dtu_qb_vs_saturn.png", "Isoform-level effect sizes: closed-form quasi-binomial vs satuRn.", "520px")}
<div><p>satuRn agrees on isoforms with ≥ 5% usage (Spearman 0.88, n = 43,479). With its regular FDR, 4,495 of 4,608 quasi-binomial calls are also significant; with its empirical-null FDR, 2,712 isoforms in 1,900 genes, all within the quasi-binomial set. For rare isoforms (&lt; ~1% usage) satuRn returns near-constant estimates with very small standard errors; these are excluded by the ΔPSI filter.</p></div></div>

<h2 id="biology">Known biology examples</h2>
<p>PSI per cell is sparse (a few reads per event), so the UMAP shows PSI smoothed over each cell's gene-expression kNN neighbourhood: (Σ inclusion reads) / (Σ event reads) over the cell and its neighbours.</p>
{img("umap_psi_known_biology.png", "kNN-smoothed PSI on the gene-expression UMAP (day 0 left cluster, day 4 right).")}
{img("psi_violin_known_biology.png", "Per-cell PSI for cells with ≥ 3 event reads; diamonds are pseudobulk PSI per day.")}
{table(kb, index=True)}
{img("psi_by_cluster_known_biology.png", "Pseudobulk PSI per Leiden cluster (≥ 20 reads), coloured by day. Tests whether shifts are uniform across subpopulations.")}
<div class="callout">
<p><b>CIRBP</b>: shifts from the proximal last exon to downstream alternative last exons and retains its 3′-proximal intron more at day 4. Temperature-dependent splicing / 3′ processing of Cirbp is documented (Gotic et al. 2016, Genes Dev), so this is a positive control for the 33.5 → 39.5 °C shift, alongside the 13-fold drop of RBM3.</p>
<p><b>SR proteins</b>: SRSF2 poison exon, SRSF5 alternative 5′ site, SRSF7 intron retention and HNRNPDL exon inclusion move together toward NMD-associated forms, while spliceosome genes go down. This fits temperature-controlled CLK/SR protein activity (Haltenhof et al. 2020, Mol Cell) and/or lower splicing demand in arrested cells. Leiden cluster 2 (a small GADD45A/PTTG1-high population, 142 day 4 and 26 day 0 cells) has the highest SRSF2 poison-exon and HNRNPDL inclusion in both days, so these events also track a stress state, not only the day.</p>
<p><b>Histones</b>: replication-dependent H1-2 and H2AC6 switch toward polyadenylated 3′ ends, expected when cells leave S phase.</p>
<p><b>Cytoskeleton</b>: TPM2 exon 6a/6b, TPM1 and CALD1 alternative first exons switch, consistent with the move from proliferating to matrix-producing cells.</p>
<p><b>RPS24</b>: inclusion of its short alternative exon rises from 23% to 96% (event PSI). The largest change in the data; whether it is driven by temperature or differentiation is not resolved.</p>
</div>
<h3>Isoform structures</h3>
<p class="muted">Top isoforms (by usage) per gene with pseudobulk usage per day; exon fill shows Δ usage (day 4 − day 0). Introns are compressed.</p>
{"".join(img(f"isoforms/{g}.png", g) for g in iso_figs)}

<h2 id="caveats">Caveats</h2>
<ul>
<li>One library per day: day is confounded with library and batch, and cells are not biological replicates. FDRs reflect cell-to-cell variability within one culture and are anti-conservative; effect sizes and biological consistency matter more.</li>
<li>Temperature and differentiation are confounded in this design; a 39.5 °C non-differentiating control would separate them.</li>
<li>Gene counts are built from pigeon-filtered isoforms, which drops ~40% of reads (mainly intra-priming).</li>
</ul>
"""

body = build()
page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>hFOB Kinnex Report</title><style>{css}</style></head><body><main>{body}</main></body></html>"""
with open(OUT, "w") as f:
    f.write(page)
print(OUT, round(os.path.getsize(OUT) / 1e6, 1), "MB")

# GitHub-rendered markdown version (GitHub shows .html files as source, not rendered)
MODE = "md"
md = build().replace("<nav>", "<p>").replace("</nav>", "</p>").replace('<div class="grid">', "<div>")
md = md.replace(' style="max-width:520px"', "").replace("</a><a", "</a> · <a").replace("</a>\n<a", "</a> · <a")
with open(f"{R}/REPORT.md", "w") as f:
    f.write(md)
print(f"{R}/REPORT.md", round(os.path.getsize(f"{R}/REPORT.md") / 1e3, 1), "kB")
