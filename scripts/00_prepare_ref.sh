#!/bin/bash
# Prepare pigeon reference: GENCODE v39 GTF, hg38 FASTA, refTSS v3.3 CAGE peaks, polyA motif list.
# usage: 00_prepare_ref.sh <ref_dir> <hg38.fa>
set -euo pipefail
R=$1; FA=$2
PB=$HOME/venv/kinnex/bin
cd "$R"
ln -sf "$FA" hg38.fa; ln -sf "$FA.fai" hg38.fa.fai
zcat gencode.v39.annotation.gtf.gz > gencode.v39.annotation.gtf
zcat refTSS_v3.3_human_coordinate.hg38.bed.gz | grep -v '^#' > refTSS_v3.3.hg38.bed  # pigeon requires the 9-column form
printf "AATAAA\nATTAAA\nAGTAAA\nTATAAA\nCATAAA\nGATAAA\nAATATA\nAATACA\nAATAGA\nAAAAAG\nACTAAA\nAAGAAA\nAATGAA\nTTTAAA\nAAAACA\nGGGGCT\n" > polyA.list.txt
$PB/pigeon prepare gencode.v39.annotation.gtf refTSS_v3.3.hg38.bed hg38.fa
