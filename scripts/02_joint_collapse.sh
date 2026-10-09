#!/bin/bash
# Merge per-sample prefixed BAMs, collapse into a joint isoform set, classify and filter with pigeon.
# usage: 02_joint_collapse.sh <work_dir> <ref_dir> <threads> <bam1> <bam2> ...
set -euo pipefail
W=$1; R=$2; T=$3; shift 3
PB=$HOME/venv/kinnex/bin
module load samtools/1.21
mkdir -p "$W"; cd "$W"
[ -s merged.bam ] || { samtools merge -f -@ $T -o merged.bam "$@"; samtools index -@ $T merged.bam; }
[ -s joint.gff ] || $PB/isoseq collapse -j $T --log-level INFO --log-file collapse.log merged.bam joint.gff
$PB/pigeon sort joint.gff -o joint.sorted.gff
$PB/pigeon classify -j $T --fl joint.abundance.txt --cage-peak $R/refTSS_v3.3.hg38.sorted.bed --poly-a $R/polyA.list.txt \
  -o joint joint.sorted.gff $R/gencode.v39.annotation.sorted.gtf $R/hg38.fa
$PB/pigeon filter -j $T joint_classification.txt --isoforms joint.sorted.gff
$PB/pigeon report joint_classification.filtered_lite_classification.txt joint_saturation.txt || true
