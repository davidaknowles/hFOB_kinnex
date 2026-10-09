#!/bin/bash
# Gene-body coverage (RSeQC) on a 2% subsample of real-cell reads per sample.
# usage: 11_genebody_coverage.sh <slurm_outputs_root> <ref_dir> <out_dir>
set -euo pipefail
SL=$1; R=$2; O=$3
module load samtools/1.21
mkdir -p $O; cd $O
bams=()
for d in day0 day4; do
  [ -s $d.sub.bam ] || samtools view -@ 4 -b -s 0.02 -d rc:1 -F 0x904 -o $d.sub.bam $SL/$d/2-SLoutputs/scisoseq.mapped.bam
  samtools index $d.sub.bam
  bams+=($O/$d.sub.bam)
done
$HOME/venv/kinnex_py/bin/geneBody_coverage.py -r $R/gencode.v39.canonical_pc.bed12 -i $(IFS=,; echo "${bams[*]}") -o genebody -f png
