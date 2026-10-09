#!/bin/bash
# Keep real-cell reads from a SMRT Link scisoseq.mapped.bam and prefix read names and CB tags with the sample label,
# so that reads/cells from different samples stay distinct after merging.
# usage: 01_prefix_cells.sh <in.mapped.bam> <sample_label> <out.bam> [threads]
set -euo pipefail
IN=$1; S=$2; OUT=$3; T=${4:-8}
module load samtools/1.21
samtools view -h -@ $T -d rc:1 "$IN" \
  | awk -v s="$S" 'BEGIN{OFS="\t"} /^@/{print; next} {$1=s"_"$1; for(i=12;i<=NF;i++) if(substr($i,1,5)=="CB:Z:"){$i="CB:Z:"s"_"substr($i,6)}; print}' \
  | samtools view -b -@ $T -o "$OUT" -
samtools index -@ $T "$OUT"
