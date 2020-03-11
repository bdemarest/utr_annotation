#! /usr/bin/bash

files=$(find ../../bam/ -name "*.bam")

for filepath in $files
do
  newfilename=$(basename $filepath _Aligned.sortedByCoord.out.bam)"_utr_test_regions.bam"
  samtools view -b -L ./utr_test_myl6_tbx5a_robo4_bed.txt $filepath > $newfilename
  samtools index $newfilename
done

