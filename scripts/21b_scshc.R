# sc-SHC (Grabski et al. 2023): significance of existing Leiden clusters (testClusters) and de novo clustering (scSHC), per day.
suppressPackageStartupMessages({library(scSHC); library(Matrix)})
args <- commandArgs(TRUE); D <- args[1]; s <- args[2]; ncores <- as.integer(args[3])
X <- as(readMM(file.path(D, paste0(s, "_counts.mtx"))), "CsparseMatrix")
rownames(X) <- readLines(file.path(D, paste0(s, "_genes.tsv")))
cells <- read.delim(file.path(D, paste0(s, "_cells.tsv")), colClasses = "character")
colnames(X) <- cells$cell
lab <- cells$leiden; names(lab) <- cells$cell
big <- names(which(table(lab) >= 20)); keep <- lab %in% big
set.seed(0)
tc <- testClusters(X[, keep], as.character(lab[keep]), alpha = 0.05, num_features = 2500, num_PCs = 30, cores = ncores)
write.table(data.frame(cell = colnames(X)[keep], leiden = lab[keep], scshc_merged = tc[[1]]),
            file.path(D, paste0(s, "_testClusters.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(D, paste0(s, "_testClusters_tree.txt"))); print(tc[[2]]); sink()
set.seed(0)
dn <- scSHC(X, alpha = 0.05, num_features = 2500, num_PCs = 30, cores = ncores)
write.table(data.frame(cell = colnames(X), leiden = lab, scshc_denovo = dn[[1]]),
            file.path(D, paste0(s, "_denovo.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(D, paste0(s, "_denovo_tree.txt"))); print(dn[[2]]); sink()
