library(readr)
library(dplyr)
require(GenomicRanges)
library(SomaticSignatures)
library(tibble)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Matching)
library(tidyr)
library(parallel)

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/PHD_PROJECT/")

# load alignability
alig = read_tsv('data/supp/wgEncodeCrgMapabilityAlign75merMask.bed', col_names = F)
colnames(alig) = c("chr", "start", "end", "alignability")

# bed files are 0-based (+1 in the start column to correct for this)
gr_alig_prev = GenomicRanges::GRanges(seqnames=alig$chr, 
              ranges = IRanges(start = alig$start+1, end = alig$end))
sum(width(gr_alig_prev))

#liftover regions
bed1 = read_tsv("../NMF2/data/FASTA_BED.ALL_GRCh37.novel_CUPs.bed", col_names = F)
bed2 = read_tsv("../NMF2/data/FASTA_BED.ALL_GRCh38.novel_CUPs.bed", col_names = F)
bed = rbind(bed1, bed2)
bed = GRanges(seqnames = bed$X1, ranges = IRanges(start = bed$X2+1 , end = bed$X3))
bed = reduce(bed)
sum(width(bed))

# remove liftover problematic regions (CUPs) from alignable regions
gr_alig_prev2 = setdiff(gr_alig_prev, bed)
sum(width(gr_alig_prev2))

# load exons
exon = read_csv('data/supp/exon_regions_processed.csv') # 0-based
exon = exon[exon$chr %in% unique(alig$chr),]
# remove exon +-2 nts
gr_exon = GenomicRanges::GRanges(seqnames = exon$chr, ranges = IRanges(start = as.numeric(exon$start)+1-2, end = as.numeric(exon$end))+2)
sum(width(gr_exon))

# remove exons from alignable regions
gr_alig = setdiff(gr_alig_prev2, gr_exon)
sum(width(gr_alig))

# load CTCFs
# fran calculated this one:
ctcf1 = read_tsv("/g/strcombio/fsupek_data/_LEGACY/epi/cohesin/CTCF_motif_matches/ren_hg19.motifs.txt")
colnames(ctcf1)[2] = "chr"
ctcf1$end = ctcf1$start + 20
ctcf1$chr = paste0("chr", ctcf1$chr)
# not sure if it is 0-based or 1-based --> I include both here to be sure
gr_ctcf1 = GenomicRanges::GRanges(seqnames = ctcf1$chr, ranges = IRanges(start = as.numeric(ctcf1$start), end = as.numeric(ctcf1$end)))
gr_ctcf1 = reduce(gr_ctcf1)

# david calculated this one with homer:
ctcf2 = read_tsv("/g/strcombio/fsupek_data/users/dmas/data/OWN_MOTIF/CTCF/20200310_CTCF_homer_test/CTCF_Zf.bed", col_names = F)
colnames(ctcf2) = c("chr", "start", "end", "type", "score", "strand")
gr_ctcf2 = GenomicRanges::GRanges(seqnames = ctcf2$chr, ranges = IRanges(start = as.numeric(ctcf2$start)+1, end = as.numeric(ctcf2$end)))
gr_ctcf2 = reduce(gr_ctcf2)

gr_ctcf_total = c(gr_ctcf1, gr_ctcf2)
gr_ctcf_total = reduce(gr_ctcf_total)
sum(width(gr_ctcf_total))

# remove CTCFs regions from alignable regions
new_gr_alig = setdiff(gr_alig, gr_ctcf_total)
sum(width(new_gr_alig))

# load ETS family
# downloaded from here: http://funseq2.gersteinlab.org/data/2.1.0
ets = read_tsv("data/supp/ETS_TFBS_ENCODE.tf.bound.union.bed", col_names = F)
colnames(ets) = c("chr", "start", "end", "name", "dot", "strand", "gene")
gr_ets = GenomicRanges::GRanges(seqnames = ets$chr, ranges = IRanges(start = as.numeric(ets$start)+1, end = as.numeric(ets$end)))
gr_ets = reduce(gr_ets)
sum(width(gr_ets))

# remove ETS TFBS regions from alignable regions
new_gr_alig2 = setdiff(new_gr_alig, gr_ets)
sum(width(new_gr_alig2))

# load APOBEC hairpins
# download from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6731024/
apo = read_tsv("/g/strcombio/fsupek_data/users/dmas/data/MISC_STRfeat/Buisson_2019_Science_TChairPins/aaw2872-Data-S2.txt")
head(apo)
apo$chr = paste0("chr", apo$chr)
apo$start = apo$pos - apo$looppos
apo$end = apo$start + apo$looplen
# not sure if it is 0-based or 1-based --> I include both here to be sure
gr_apo = GenomicRanges::GRanges(seqnames = apo$chr, ranges = IRanges(start = as.numeric(apo$start), 
                                                                     end = as.numeric(apo$end)))
gr_apo = reduce(gr_apo)
sum(width(gr_apo))

# remove APOBEC hairpin regions from alignable regions
new_gr_alig3 = setdiff(new_gr_alig2, gr_apo)
sum(width(new_gr_alig3))

# black list regions!!
bl = read_tsv("../NMF2/data/hg19-blacklist.v2.bed.gz", col_names = F)
head(bl)
table(bl$X4)
colnames(bl) = c("chr", "start", "end", "sig")
gr_bl= GenomicRanges::GRanges(seqnames = bl$chr, 
       ranges = IRanges(start = as.numeric(bl$start)+1, 
                        end = as.numeric(bl$end)))
gr_bl = reduce(gr_bl)
sum(width(gr_bl))

# remove APOBEC hairpin regions from alignable regions
new_gr_alig4 = setdiff(new_gr_alig3, gr_bl)
sum(width(new_gr_alig4))

# write dataframe
head(new_gr_alig4)
df = data.frame(chr = seqnames(new_gr_alig4), start = start(new_gr_alig4), 
                end = end(new_gr_alig4))
sum(df$end - df$start +1)
write.csv(df, "data_v2/filters/processed_alignability_LIFTOVERCUPs_BLACKLIST_crg75_without_exon_ctcf_ets_apobec_hairpins.csv", row.names = F)
