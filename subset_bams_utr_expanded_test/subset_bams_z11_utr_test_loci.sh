#! /usr/bin/bash

files=$(find ../../bam_z11/ -name "*.bam")

for filepath in $files
do
  newfilename=$(basename $filepath _Aligned.sortedByCoord.out.bam)"_utr_test_regions.bam"
  samtools view -b -L ./20200311_testgenesfor3UTR_expand10000bp.bed $filepath > $newfilename
  samtools index $newfilename
done

