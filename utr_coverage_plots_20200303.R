# Look at 3' UTR coverage for small genes-of-interest list.
# Evaluate heuristic 3' UTR extension.
# Test using z10 STAR bam files.


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

gtf[, gene_biotype:=gsub('^.+; gene_biotype \\"(.+?)\\"; .+$', "\\1", attributes)]

gtf[, transcript_biotype:=gsub('^.+; transcript_biotype \\"(.+?)\\"; .+$', "\\1", attributes)]
gtf[nchar(transcript_biotype) > 34, transcript_biotype:=NA_character_]


addmargins(table(nchar(gtf$transcript_id), useNA="always"))
# 18    <NA>     Sum 
# 1129345   32520 1161865 

dim(gtf)
# [1] 1161865      12

# Table of annotated stop_codon positions from 
subtab = gtf[feature == "stop_codon",
             list(sequence, start, end, strand, gene_id,
                  transcript_id, gene_name, gene_biotype, transcript_biotype)]

# How many gene_name have multiple gene_id?
gid_tab = gtf[feature == "gene", list(n_gene_id=length(unique(gene_id))), by=gene_name]


#===============================================================================






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


