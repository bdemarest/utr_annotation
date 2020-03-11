#! /bin/bash

# Compute coverage of regions of bam file using bedtools.

files=$(find . -name "*.bam")

for bamfile in $files
    do
        newfilename_pos=$(basename $filepath .bam)"_pos_coverage.txt.gz"
        newfilename_neg=$(basename $filepath .bam)"_neg_coverage.txt.gz"

        bedtools coverage \
                 -a utr_test_regions_bed_20200311.txt \
                 -b $bamfile \
                 -d -bed -s | gzip > $newfilename_pos

        bedtools coverage \
                 -a utr_test_regions_bed_20200311.txt \
                 -b $bamfile \
                 -d -bed -S | gzip > $newfilename_neg
    done
