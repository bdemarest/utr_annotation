# Look at 3' UTR coverage for small genes-of-interest list.
# Evaluate heuristic 3' UTR extension.
# Test using z11 STAR bam files.

# Github notes.
# When you have an existing Rstudio project that you want to bring 
# into git and github:
# (Assuming you have GITHUB_PAT in your ~/.Renviron folder)

# (1) Add project to local git (Git tab will now appear in Rstudio window).
# usethis::use_git()
# (2) Create github repository, commit, and push to github.
# usethis::use_github()


library(data.table)
library(ggplot2)
library(biomaRt)
library(Rsamtools)
library(fs)
library(here)



# Load bed file of genes of interest.
goi = fread(here("subset_bams_utr_expanded_test", 
                 "20200311_testgenesfor3UTR_realpositions.bed"))

setnames(goi, c("chromosome", "start", "end", "gene_id_combined"))

goi[, c("external_gene_name", "ensembl_gene_id"):=tstrsplit(gene_id_combined, "_", fixed=TRUE)]

#===============================================================================
# Can we get STOP codon locations from gtf file?

if (file_exists(here("gtf_data_z11_r99.txt.gz"))) {
  # Load pre-parsed gtf data (tab-delimited format) from gzipped file.
  gtf = fread(here("gtf_data_z11_r99.txt.gz"))
} else {
  # Or, load from original gtf file, do all parsing steps, and save to .txt.gz.
  gtf = fread(here::here("Danio_rerio.GRCz11.99.gtf"), header=FALSE)
  
  gtf_col_names = c(
    "sequence",
    "source",
    "feature", 
    "start", 
    "end",
    "score", 	
    "strand",	
    "phase",
    "attributes"
  )
  
  setnames(gtf, gtf_col_names)
  
  gtf[, gene_id:=gsub("^.+(ENSDARG\\d{11}).+$", "\\1", attributes)]
  gtf[, transcript_id:=gsub("^.+(ENSDART\\d{11}).+$", "\\1", attributes)]
  gtf[nchar(transcript_id) > 18, transcript_id:=NA_character_]
  gtf[, gene_name:=gsub('^.+; gene_name \\"(.+?)\\"; .+$', "\\1", attributes)]
  gtf[, gene_biotype:=gsub('^.+; gene_biotype \\"(.+?)\\";.*$', "\\1", attributes)]
  gtf[, transcript_biotype:=gsub('^.+; transcript_biotype \\"(.+?)\\"; .+$', "\\1", attributes)]
  gtf[nchar(transcript_biotype) > 34, transcript_biotype:=NA_character_]
  
  addmargins(table(nchar(gtf$transcript_id), useNA="always"))
  # 18    <NA>     Sum 
  # 1129345   32520 1161865 
  
  dim(gtf)
  # [1] 1161865      12
  
  fwrite(gtf, file=here("gtf_data_z11_r99.txt.gz"), sep="\t")
}

format(object.size(x=gtf) ,units="MB")

#-------------------------------------------------------------------------------
# Table of annotated stop_codon positions from 
stop_tab = gtf[feature == "stop_codon",
               list(sequence, stop_codon_start=start, stop_codon_end=end, 
                    strand, gene_id, transcript_id, gene_name, 
                    gene_biotype, transcript_biotype)]

utr_tab = gtf[feature == "three_prime_utr",
              list(sequence, three_prime_utr_start=start, three_prime_utr_end=end, 
                   strand, gene_id, transcript_id, gene_name, 
                   gene_biotype, transcript_biotype)]

subtab = merge(x=stop_tab, y=utr_tab, all.x=TRUE, all.y=TRUE)

#===============================================================================
# Another way to load gtf info. Loads as GRanges oject.

library(rtracklayer)

# Also works directly on .gz file.
gtf_gr = import(here("Danio_rerio.GRCz11.99.gtf"))

format(object.size(x=gtf_gr) ,units="MB")

#===============================================================================
#----------------
# n = 24007 unique ensembl gene ids have associated stop codon.
length(unique(stop_tab$gene_id))
# Out of 25432 protein coding genes in the full gtf file.
length(unique(gtf[feature == "gene" & gene_biotype == "protein_coding", gene_id]))


# We only have 39559 transcripts with associated stop codons.
length(unique(stop_tab$transcript_id))
# [1] 39559
length(unique(gtf$transcript_id))
# [1] 59877

# Example: ptpn12 ENSDARG00000102141 is annotated as protein coding,
# and has a refseq predicted protein sequence, but has no
# ensembl annotated utrs, nor stop codon.

# missing transcript ids.
mtxid_vec = setdiff(unique(gtf$transcript_id), unique(stop_tab$transcript_id))

mtx_tab = gtf[transcript_id %in% mtxid_vec]
# ?? Why are rows with NA transcript id showing up in this table??

#----------------

# > dim(subtab)
# [1] 43636    11
# > length(unique(subtab$transcript_id))
# [1] 39689
# > dim(stop_tab)
# [1] 39635     9
# > length(unique(stop_tab$transcript_id))
# [1] 39559
# > dim(utr_tab)
# [1] 34772     9
# > length(unique(utr_tab$transcript_id))
# [1] 30928
# > 43636 - 39689
# [1] 3947

# Many transcript ids have more than one associated utr.
# A few transcript ids have more than one associated stop codon.

#----------------

# How many gene_name have multiple gene_id?
gid_tab = gtf[feature == "gene", list(n_gene_id=length(unique(gene_id))), by=gene_name]
table(gid_tab$n_gene_id)
#     1     2     3     4     5     6     8     9    10    18    19    24 
# 31256   474    47     8     5     2     2     1     2     1     1     1 

# For example, 
gid_tab[n_gene_id == 24]
#    gene_name n_gene_id
# 1:  hist1h4l        24
#===============================================================================
# Move forward with subtab table of stop codon and utr positions.


goi_tab = subtab[gene_name %in% goi$external_gene_name]

# Add column of unique row index values. (Because tx id is not a unique id.)
goi_tab[, row_index:=.I]


#-------------------------------------------------------------------------------
# Retrieve coverage from bam files using Rsamtools.

# Some good ideas here:
# https://support.bioconductor.org/p/69620/
# Also, BamViews strategy from 
# http://bioconductor.org/help/course-materials/2010/BioC2010/Exercises-Rsamtools.pdf
# Another option:
# https://kasperdanielhansen.github.io/genbioconductor/html/Rsamtools.html

bam_file_names = list.files(path=here("subset_bams_utr_expanded_test"), 
                            pattern="*.bam$", full.names=TRUE, recursive=TRUE)


gnomex_id_vec = gsub("^.+\\D(\\d{1,5}X\\d{1,3}).+$", "\\1", bam_file_names)

# Make GRanges object defining stop codon positions for all genes-of-interest.
# 
goi_ranges = GRanges(seqnames=goi_tab$sequence,
                     strand=goi_tab$strand,
                     ranges=IRanges(start=goi_tab$stop_codon_start,
                                    end=goi_tab$stop_codon_start),
                     ensembl_gene_id=goi_tab$gene_id,
                     ensembl_transcript_id=goi_tab$transcript_id,
                     external_gene_name=goi_tab$gene_name)
                     
# Add 10,000 bases 3' of stop codon using start=FALSE argument.
goi_ranges = flank(goi_ranges, start=FALSE, width=10000)

# https://bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
# One way of understanding a BamViews instance is as a rectangular data
# structure. The columns represent BAM files (e.g., distinct samples). The rows
# represent ranges (i.e., genomic coordinates).

# Check available methods.
methods(class="BamViews")

bv = BamViews(bam_file_names, bamSamples=DataFrame(gnomex_id=gnomex_id_vec),
              bamRanges=goi_ranges)
names(bv) = gnomex_id_vec

params = ScanBamParam(which=goi_ranges, what=c("rname", "qname", "pos", "qwidth"))

# res is a nested list containing info from each bam file,
# for each region of interest, on each read select by the GRanges criteria.
res = scanBam(bv, param=params)

# Need to convert read positions to coverage
# ?? coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))





#------------------------------- 
# Another possible method using Rsamtools::applyPileups
# Get the base counts.
goi_param = ApplyPileupsParam(which=goi_ranges, what="seq",
                              minMapQuality=13, maxDepth=1000L)

res_list = list()
for (i in seq_along(bam_files)) {
  bamfile = bam_files[i]
  gnomex_id = gnomex_id_vec[i]
  bam_id = bam_id_vec[i]
  tmp_pileup = applyPileups(files=PileupFiles(bamfile),
                            FUN=function(x) x,
                            param=chr3_param)
  
  tmp = data.table(t(tmp_pileup[[1]][["seq"]][1:4, 1, ]))
  set(tmp, j="pos", value=tmp_pileup[[1]][["pos"]])
  set(tmp, j="gnomex_id", value=gnomex_id)
  set(tmp, j="bam_id", value=bam_id)
  res_list[[i]] = tmp
}

full_tab = rbindlist(res_list)





# https://homolog.us/Bioconductor/Rsamtools.html
## Any combination of the filtering criteria is possible: let's say we
## want a "coverage pileup" that only counts reads with mapping
## quality >= 13, and bases with quality >= 10 that appear at least 4
## times at each genomic position
p_param <- PileupParam(distinguish_nucleotides=FALSE,
                       distinguish_strands=FALSE,
                       min_mapq=13,
                       min_base_quality=10,
                       min_nucleotide_depth=4)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
head(res)
res <- res[, c("pos", "count")] ## drop seqnames and which_label cols
plot(count ~ pos, res, pch=".")








# Detecting 'splice-compatible' alignments.
# http://bioconductor.org/packages/release/bioc/vignettes/GenomicAlignments/inst/doc/OverlapEncodings.pdf






# List of bam, bai, and gtf files some of which are required
# but not included in github.

list.files(recursive=TRUE)
# [1] "20200303_testgenesfor3UTR.txt"                                            
# [2] "Danio_rerio.GRCz11.99.gtf"                                                
# [3] "gtf_data_z11_r99.txt.gz"                                                  
# [4] "subset_bams_utr_expanded_test/14893X1_utr_test_regions.bam"               
# [5] "subset_bams_utr_expanded_test/14893X1_utr_test_regions.bam.bai"           
# [6] "subset_bams_utr_expanded_test/14893X2_utr_test_regions.bam"               
# [7] "subset_bams_utr_expanded_test/14893X2_utr_test_regions.bam.bai"           
# [8] "subset_bams_utr_expanded_test/14893X3_utr_test_regions.bam"               
# [9] "subset_bams_utr_expanded_test/14893X3_utr_test_regions.bam.bai"           
# [10] "subset_bams_utr_expanded_test/14893X4_utr_test_regions.bam"               
# [11] "subset_bams_utr_expanded_test/14893X4_utr_test_regions.bam.bai"           
# [12] "subset_bams_utr_expanded_test/14893X5_utr_test_regions.bam"               
# [13] "subset_bams_utr_expanded_test/14893X5_utr_test_regions.bam.bai"           
# [14] "subset_bams_utr_expanded_test/20200311_testgenesfor3UTR_expand10000bp.bed"
# [15] "subset_bams_utr_expanded_test/20200311_testgenesfor3UTR_realpositions.bed"
# [16] "subset_bams_utr_expanded_test/subset_bams_z11_utr_test_loci.sh"           
# [17] "subset_bams_utr_test/14893X1_utr_test_regions_neg_coverage.txt.gz"        
# [18] "subset_bams_utr_test/14893X1_utr_test_regions_pos_coverage.txt.gz"        
# [19] "subset_bams_utr_test/14893X1_utr_test_regions.bam"                        
# [20] "subset_bams_utr_test/14893X1_utr_test_regions.bam.bai"                    
# [21] "subset_bams_utr_test/14893X2_utr_test_regions_neg_coverage.txt.gz"        
# [22] "subset_bams_utr_test/14893X2_utr_test_regions_pos_coverage.txt.gz"        
# [23] "subset_bams_utr_test/14893X2_utr_test_regions.bam"                        
# [24] "subset_bams_utr_test/14893X2_utr_test_regions.bam.bai"                    
# [25] "subset_bams_utr_test/compute_bedtools_coverage_test_loci.R"               
# [26] "subset_bams_utr_test/compute_bedtools_coverage_test_loci.sh"              
# [27] "subset_bams_utr_test/subset_bams_utr_test_loci.sh"                        
# [28] "subset_bams_utr_test/utr_test_myl6_tbx5a_robo4_bed.txt"                   
# [29] "subset_bams_utr_test/X1_pos_coverage_test_regions_20200311.pdf"           
# [30] "subset_bams_utr_test/X1_pos_coverage_test_regions_histograms_20200311.pdf"
# [31] "utr_annotation.Rproj"                                                     
# [32] "utr_coverage_plots_20200303.R"                                            
# [33] "utr_coverage_plots_20200325.R"                                            
# [34] "utr_test_regions_bed_20200311.txt"



