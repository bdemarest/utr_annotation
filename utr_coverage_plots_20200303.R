# Look at 3' UTR coverage for small genes-of-interest list.
# Evaluate heuristic 3' UTR extension.
# Test using z10 STAR bam files.

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
library(here)
library(Gviz)
library(GenomicAlignments)
library(rtracklayer)
library(biomaRt)


library(Rsamtools)


bam_file_names = list.files(path=here("subset_bams_utr_test"), 
                            pattern="*.bam$", full.names=TRUE, recursive=TRUE)

# Load bed file of genes of interest.
goi = fread(here("subset_bams_utr_test", "utr_test_myl6_tbx5a_robo4_bed.txt"))

setnames(goi, c("chromosome", "start", "end", "external_gene_name"))


# Load biomart annotation.
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
               dataset="drerio_gene_ensembl",
               host="jan2020.archive.ensembl.org")


mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
               dataset="drerio_gene_ensembl",
               host="jan2020.archive.ensembl.org",
               mirror="useast")





annot = as.data.table(getBM(mart=mart,
                               attributes=c("ensembl_gene_id",
                                            "ensembl_transcript_id",
                                            "external_gene_name",
                                            "chromosome_name",
                                            "start_position",
                                            "end_position",
                                            "transcript_start",
                                            "transcript_end",
                                            "3_utr_start",
                                            "3_utr_end",
                                            "gene_biotype",
                                            "transcript_biotype")))

# Build a new table
# (1) If a tx_id has only one row, keep it.
# (2) If a tx_id has two or more rows, throw out any that have NA values for utr start.
# (3) If a tx_id still has more than one row, randomly select one of the rows.

annot_1 = annot[is.na(`3_utr_start`)]
annot_1 = annot_1[!duplicated(annot_1, )]



tx_count_tab = annot[, list(txid_row_count =length(external_gene_name)),
                            
                     by=ensembl_transcript_id]






#-------------------------------------------------------------------------------
# Retrieve coverage from bam files using Rsamtools.


gnomex_id_vec = gsub("^.+\\D(\\d{1,5}X\\d{1,3}).+$", "\\1", bam_files)
bam_id_vec = gsub("Aligned.sortedByCoord.out.bam", "", basename(bam_files), fixed=TRUE)

chr3_len = 63268876

goi_ranges = GRanges(seqnames="3",
                     ranges=IRanges(start=1, end=chr3_len))

goi_param = ApplyPileupsParam(which=goi_ranges, what="seq",
                               minMapQuality=13, maxDepth=1000L)

# Get the base counts.

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




#===============================================================================
# Can we get STOP codon locations from gtf file?

library(data.table)
library(here)

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

# Load bed file of test genes of interest.
# (These are the only regions with reads in the subsetted test bam files).
goi = fread(here::here("subset_bams_utr_test", "utr_test_myl6_tbx5a_robo4_bed.txt"))
setnames(goi, c("chromosome", "start", "end", "external_gene_name"))

goi_tab = subtab[gene_name %in% goi$external_gene_name]
# Add column of unique row index values.
goi_tab[, row_index:=.I]

# Try bedtools coverage.
# Make a regions bed file.

# Bed file header:
# chromosome, start, end, name, score, strand

# Only used for pulling coverage out of bam files using bedtools.
# (Can also use this exact bed file to subset the bams using samtools).
bed_tab = goi_tab[, list(chromosome=sequence,
                         start=min(stop_codon_start, three_prime_utr_start) - 1e4,
                         end=max(stop_codon_end, three_prime_utr_end) + 1e4,
                         name=paste(gene_name, transcript_id, sep="_"),
                         score=0L,
                         strand=strand),
                  by=row_index]

fwrite(bed_tab[, -"row_index"], file=here::here("utr_test_regions_bed_20200311.txt"), 
       sep="\t", col.names=FALSE)


#-------------------------------------------------------------------------------
# Load coverage data (from bedtools coverage).
library(data.table)
library(here)
library(ggplot2)

tmp = fread(here::here("subset_bams_utr_test", 
                       "14893X1_utr_test_regions_pos_coverage.txt.gz"))


p1 = ggplot(tmp, aes(x=V7, y=V8)) +
     geom_line() +
     facet_wrap(~ V4, scales="free", ncol=2)

ggsave(file=here("subset_bams_utr_test",
                 "X1_pos_coverage_test_regions_20200311.pdf"),
       plot=p1, width=10, height=12)


h1 = ggplot(tmp, aes(x=V8)) +
  geom_histogram(colour="grey30", fill="grey80", bins=100) +
  facet_wrap(~ V4, scales="free", ncol=2)

ggsave(file=here("subset_bams_utr_test",
                 "X1_pos_coverage_test_regions_histograms_20200311.pdf"),
       plot=h1, width=10, height=12)















#-------------------------------------------------------------------------------
# From https://bioinformatics.stackexchange.com/questions/525/bam-to-bigwig-without-intermediary-bedgraph

## read in BAM file (use readGAlignmentPairs for paired-end files)
gr <- readGAlignments(bam_file_n)

## convert to coverages
gr.cov <- coverage(gr)

## export as bigWig
export.bw(gr.cov,'sample1.bigwig')

#-------------------------------------------------------------------------------
# From https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html


dtrack_X1 = DataTrack(range=bam_file_names[1], type="l", name="Coverage")

class(dtrack_X1)

dtrack_X1



plotTracks(dtrack_X1)


plotTracks(dtrack_X1, chromosome=goi[1, chromosome], from=goi[1, start], to=goi[1, end])


bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
                                    name="ENSEMBL", biomart=bm)
plotTracks(biomTrack)


