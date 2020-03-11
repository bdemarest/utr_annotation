# R script to run bedtools commands. (something wrong with bash script).

library(data.table)
library(parallel)

tab = data.table(bamfile_path=list.files(pattern=".bam$", path="."))
tab[, outfile_name_pos:=paste(gsub(pattern="\\.bam$",
                                   replacement="",
                                   basename(bamfile_path)), 
                              "_pos_coverage.txt.gz", sep="")]


tab[, outfile_name_neg:=paste(gsub(pattern="\\.bam$",
                                   replacement="",
                                   basename(bamfile_path)), 
                              "_neg_coverage.txt.gz", sep="")]

command_vec = paste(
    "bedtools coverage",
    "    -a utr_test_regions_bed_20200311.txt",
    "    -b", tab$bamfile_path,
    "    -d",
    "    -bed",
    "    -s",
    "| gzip >", tab$outfile_name_pos, "\n",
     
    "bedtools coverage",
    "    -a utr_test_regions_bed_20200311.txt",
    "    -b", tab$bamfile_path,
    "    -d",
    "    -bed",
    "    -S",
    "| gzip >", tab$outfile_name_neg
)

cl = makeCluster(8)

parLapply(cl=cl, X=command_vec, fun=system)

stopCluster(cl)




