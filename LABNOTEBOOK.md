# hFOB Kinnex single-cell long-read lab notebook

## Data

PacBio Kinnex single-cell (10x Genomics 3' v3.x, 12 bp UMI), Revio, one library per time point of hFOB 1.19 osteoblast differentiation: day 0 and day 4. Inputs are SMRT Link outputs (skera → isoseq tag/refine/correct/groupdedup → pbmm2 to hg38, GENCODE v39 → isoseq collapse → pigeon) run separately per sample. hFOB 1.19 carries a temperature-sensitive SV40 large T antigen (tsA58): cells proliferate at 33.5 °C and differentiate at 39.5 °C, where large T is inactive.

## Pipeline

1. `scripts/00_prepare_ref.sh` pigeon reference: hg38, GENCODE v39 GTF, refTSS v3.3 (CAGE) peaks (pigeon requires the original 9-column BED), SQANTI polyA motif list.
2. `scripts/01_prefix_cells.sh` keep real-cell reads (`rc:1`) from each per-sample mapped BAM and prefix read names and `CB` tags with the sample label, so that cells and reads from different samples stay distinct after merging.
3. `scripts/02_joint_collapse.sh` merge, `isoseq collapse` (defaults) on the merged BAM to obtain one joint isoform set, `pigeon classify` (CAGE, polyA motifs) and `pigeon filter` (intra-priming, RT switching, junction coverage). Per-sample SMRT Link runs give non-matching PB IDs, so the joint collapse is needed for any isoform comparison between days.
4. `scripts/06_build_counts.py` (`hfob_kinnex/counts.py`) cell × isoform counts from read → cell (BAM `CB`) and read → isoform (collapse group file) for pigeon-filtered isoforms; gene counts are sums over isoforms of each `associated_gene`. Each read is already one UMI-deduplicated molecule.
5. QC: `03_read_qc.py`/`12_read_qc_plots.py` (read lengths), `04_knee_plot.py`, `11_genebody_coverage.sh`/`13_genebody_plot.py` (RSeQC gene-body coverage on 2% of reads over GENCODE Ensembl_canonical protein-coding transcripts), `14_read_categories.py` (per-read pigeon category and filter reason).
6. `07_sc_analysis.py` (`hfob_kinnex/sc.py`) per-sample MAD cell filters (3 MAD, low UMIs/genes, high mito); scanpy normalisation, seurat_v3 HVG, PCA, UMAP, Leiden in three feature spaces: gene expression, isoform expression (isoforms detected in ≥ 20 cells) and isoform usage. Usage features: for the 3000 most expressed genes detected in ≥ 20% of cells with ≥ 2 isoforms of pooled usage ≥ 5%, $u_{ci} = x_{ci} / \sum_{j \in g(i)} x_{cj} - \bar u_i$ where $x_{ci}$ is the count of isoform $i$ in cell $c$, $g(i)$ the isoforms of its gene and $\bar u_i$ the pooled usage; cells with zero gene counts get $u_{ci}=0$.
7. `08_events.py` (`hfob_kinnex/events.py`, `hfob_kinnex/stats.py`) alternative events on isoforms with ≥ 10 reads:
   - SUPPA2 local events (SE, A5, A3, MX, RI, AF, AL; strict boundaries) from the expressed-isoform GTF.
   - TSS and 3' end clusters per gene (single linkage, 50 nt for TSS, 100 nt for 3' ends); a gene has alt TSS / alt polyA if ≥ 2 clusters each have pooled usage ≥ 10% and ≥ 20 reads. Clusters are classified as distinct terminal exons (non-overlapping first/last exons) or tandem (same exon). "Supported" clusters have ≥ 50% of reads from isoforms whose TSS lies in a refTSS peak (TSS) or whose 3' end has an upstream polyA signal (3' end).
   - Detection: event with ≥ 20 pooled reads and 0.05 ≤ PSI ≤ 0.95 (0.1–0.9 for end clusters).
8. Differential usage, day4 vs day0, cells as observations. For each event / cluster / isoform, $y_c$ = inclusion (or isoform) reads and $n_c$ = event (or gene) total reads in cell $c$, $\mathrm{logit}\, p_c = a + b\, d_c$ with $d_c \in \{0,1\}$ the day indicator, and $\mathrm{Var}(y_c) = \phi n_c p_c (1-p_c)$. This is a quasi-binomial GLM with one binary covariate; the MLE is the pooled per-day proportion, $\phi$ is the Pearson $\chi^2$ over cells with $n_c>0$ divided by the residual df (floored at 1), and $b$ is tested by Wald. Implemented in closed form, vectorised over features (`quasibinomial_two_group`). Significant if BH FDR < 0.05 and |ΔPSI| ≥ 0.1. Requires ≥ 20 cells and ≥ 50 reads per day. satuRn (quasi-binomial GLM with empirical Bayes dispersion, `10a_export_dtu.py`, `10b_saturn.R`) is run on the same isoform set as a check.
9. `09_gene_de.py` cell-level Wilcoxon (scanpy) with pseudobulk log2 CPM fold change; significant if FDR < 0.05, |log2FC| ≥ 0.5 and detected in ≥ 10% of cells in either day. GSEA prerank (gseapy) on pseudobulk log2FC against MSigDB Hallmark 2020 and GO BP 2023.
10. `15_isoform_examples.py` (`hfob_kinnex/plots.py`) isoform structure plots with per-day usage; `16_cluster_markers.py` within-day Leiden markers.

Caveat for all tests: there is one library per day, so day is confounded with library/batch, and cells are not biological replicates. p-values reflect cell-to-cell variability within one culture and are anti-conservative for inference about the differentiation process. Effect sizes (log2FC, ΔPSI) and biological consistency carry more weight than FDR.

## Findings (2026-10-09)

### Sequencing and read QC
- SMRT Link: 10,258 (day 0) and 5,939 (day 4) cells; 140M and 131M S-reads; ~16-mer arrays (HiFi array reads ~11.5–12.3 kb). Knee plots show clean cell / background separation.
- S-reads median 642 / 692 nt; mapped real-cell reads median 515 / 556 nt (10th–90th pct ~300–970 nt). This is short for full-length 10x cDNA; ~45% of reads are mono-exonic.
- Joint collapse: 196.6M real-cell reads in 16,180 cells, 7.67M raw isoforms; 1.28M pass pigeon filter (128k FSM, 265k ISM, 152k NIC, 576k NNC). By read (QC-pass cells): 60% pass, 19–20% intra-priming, 17% low coverage/non-canonical junction, 3% RT switching. Of passing reads, 66% are full-splice matches.
- 5'/3' bias: gene-body coverage on canonical protein-coding transcripts peaks at the 5' end and declines toward the 3' end (10th/90th percentile coverage ratio 1.36). Among passing reads assigned to known transcripts, 5'-fragments (ISMs missing 3' exons) are 8.5% / 7.3% and 3'-fragments (missing 5' exons) are 3.0%. So molecules mostly reach the TSS; the 3' deficit comes from early 3' ends (internal oligo-dT priming not caught by the filter, and proximal polyA sites upstream of the long annotated 3' UTR). Both days have the same profile.

### Cells
- QC-pass: 10,065 day 0 and 5,653 day 4 cells. Median UMIs 5,106 / 6,815, genes 1,843 / 2,257, isoforms 2,301 / 2,887 per cell (after pigeon filtering; SMRT Link medians before filtering are 8,718 / 12,589). Mitochondrial fraction is low (< 1%).
- Days separate completely in gene expression, isoform expression and isoform usage UMAPs. Within days, structure is mostly cell cycle (day 0: 62% G2M, 31% S; day 4: 46% G1, 12% G2M), plus a small interferon-response cluster in day 0 (ISG15, IFIT1–3; 213 cells) and a POSTN/FN1/COL4A1/INHBA matrix-producing subset at day 4 with the highest osteoblast score.

### Differential expression (day 4 vs day 0)
1,501 up and 1,553 down genes. The changes match the tsA58 temperature shift and early osteoblast differentiation:
- p53 reactivation after large T inactivation: Hallmark p53 pathway is the top up-regulated set (NES 2.6); CDKN1A (+3.6 log2FC, detected in 99% of day 4 cells), MDM2 (+3.1), GDF15 (+5.3), TP53I3 (+3.6), BAX (+1.9).
- Proliferation arrest: G2M checkpoint, E2F and MYC targets are the top down-regulated sets; CCNB1 −2.9, TOP2A −1.8, MKI67 −1.5; replication-dependent histones strongly down.
- Osteoblast program: ALPL +3.5 (detected in 0.8% → 17% of cells), SPP1 +3.6, BGLAP +1.9 (low), COL1A2 +2.1, COL1A1 +0.8, POSTN +2.0, RUNX2 +0.7; EMT/ECM organisation up.
- Temperature: cold-inducible RBM3 −3.7 log2FC and CIRBP −0.5, consistent with 33.5 → 39.5 °C. Classical heat-shock genes are not induced (HSPA1A/B, DNAJB1 unchanged), so cells have adapted by day 4.
- Inflammatory / SASP-like sets up (SAA1, CXCL2, IL-6/JAK/STAT3, TNF/NF-κB), consistent with p53-driven arrest. GO mRNA splicing / spliceosome genes are down (NES −2.6).

### Alternative isoform events
Detectable events (≥ 20 pooled reads, 0.05 ≤ PSI ≤ 0.95): SE 11,701 (5,126 genes), A5 7,939, A3 6,746, MX 3,145, RI 3,030, AF 53,064 (3,754 genes), AL 31,314 (5,004 genes). Genes with ≥ 2 used TSS clusters: 8,185 (895 when every cluster must be CAGE-supported); with ≥ 2 used 3' end clusters: 9,509 (3,314 with polyA-signal support). The large drop under CAGE support indicates many apparent alternative TSSs are 5'-truncated molecules; the supported numbers are the conservative estimate. SUPPA AF/AL counts are inflated by combinatorial pairs of isoforms and the cluster-based counts are preferred.

Differential (FDR < 0.05, |ΔPSI| ≥ 0.1): SE 1,845 events / 1,258 genes; A5 954 / 647; A3 681 / 516; MX 491 / 309; RI 440 / 300; TSS clusters 3,022 / 1,823 genes; 3' end clusters 3,879 / 2,470 genes; isoform usage 4,608 isoforms / 3,215 genes.

Biological checks:
- CIRBP shifts from its proximal last exon (71% of reads at day 0 → 38% at day 4) to isoforms splicing into downstream alternative last exons (9% → 43%), and its 3'-proximal intron is retained more (PSI 0.59 → 0.90). Temperature-dependent splicing / 3' end processing of Cirbp is documented (Gotic et al. 2016 Genes Dev); this is a positive control for the temperature shift.
- SR protein autoregulatory exons shift together: SRSF2 poison exon, SRSF3 exon 4, SRSF5 A5, SRSF7 RI, SRSF10 A3 and HNRNPDL (ΔPSI +0.29) all move towards NMD-associated forms at day 4, alongside lower expression of splicing genes. This fits temperature-dependent CLK/SR protein activity (Haltenhof et al. 2020 Mol Cell) and/or reduced splicing demand in arrested cells.
- Replication-dependent histones (H1-2, H2AC6) switch to the longer, polyadenylated 3' isoform (H1-2: 25% → 66%), as expected when cells leave S phase and stem-loop processing declines.
- Cytoskeletal genes TPM1 (alternative first exons), TPM2 (MX exons 6a/6b, 46% → 9%) and CALD1 (alternative first exons) change usage, consistent with the shift from proliferating to matrix-producing, contractile-like cells.
- RPS24 switches to the isoform including its short alternative exon (21% → 91% of reads). This is the largest change in the data. RPS24 alternative exons are known to be tissue-regulated; the cause here (temperature vs differentiation) is not resolved.

## Open questions
- Differentiation protocol: whether day 4 was at 39.5 °C only or also in osteogenic medium. Data are consistent with the temperature shift; the osteoblast changes are modest (ALPL in 17% of cells) as expected at day 4.
- Separating temperature from differentiation effects would need a 39.5 °C non-differentiating control or a parental line without tsA58.
- Replicate libraries per time point would allow pseudobulk tests with proper biological variance.
