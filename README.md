# hFOB Kinnex single-cell long-read analysis

Processing and analysis of PacBio Kinnex single-cell (10x 3') long-read data from hFOB 1.19 osteoblasts at day 0 and day 4 of differentiation: joint isoform discovery across samples, read and cell QC, embeddings on gene expression, isoform expression and isoform usage, detection of alternative splicing / TSS / polyA events, and day 4 vs day 0 differential expression and isoform usage.

- `hfob_kinnex/` reusable functions (counts, single-cell helpers, event definitions, quasi-binomial usage test, plots)
- `scripts/` numbered pipeline steps, run in order
- `results/figures`, `results/tables` outputs
- `LABNOTEBOOK.md` methods and findings

Tools: isoseq (collapse), pigeon (classify/filter), samtools, SUPPA2, RSeQC, scanpy, gseapy, satuRn. Python dependencies are installed with uv; isoseq and pigeon come from bioconda.
