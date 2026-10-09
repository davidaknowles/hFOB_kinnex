# satuRn differential transcript usage day4 vs day0, cells as observations.
suppressPackageStartupMessages({library(satuRn); library(SummarizedExperiment); library(Matrix); library(BiocParallel)})
args <- commandArgs(TRUE); D <- args[1]; ncores <- as.integer(args[2])
counts <- as(readMM(file.path(D, "counts.mtx")), "CsparseMatrix")
tx <- read.delim(file.path(D, "txinfo.tsv")); cells <- read.delim(file.path(D, "cells.tsv"))
rownames(counts) <- tx$isoform; colnames(counts) <- cells$cell
txInfo <- data.frame(isoform_id = tx$isoform, gene_id = tx$gene, row.names = tx$isoform)
cd <- data.frame(group = factor(cells$sample, levels = c("day0", "day4")), row.names = cells$cell)
se <- SummarizedExperiment(assays = list(counts = as.matrix(counts)), colData = cd, rowData = txInfo)
se <- fitDTU(se, formula = ~ 0 + group, parallel = TRUE, BPPARAM = MulticoreParam(ncores), verbose = TRUE)
design <- model.matrix(~ 0 + group, data = cd)
L <- matrix(c(-1, 1), ncol = 1, dimnames = list(colnames(design), "day4_vs_day0"))
se <- testDTU(se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)
res <- rowData(se)[["fitDTUResult_day4_vs_day0"]]
res <- cbind(isoform = rownames(se), gene = rowData(se)$gene_id, as.data.frame(res))
write.table(res, gzfile(file.path(D, "saturn_day4_vs_day0.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
