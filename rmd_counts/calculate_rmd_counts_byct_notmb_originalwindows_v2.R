#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(require(GenomicRanges))
library(parallel)

setwd('/g/strcombio/fsupek_home/msalvadores/Documents/RMDSigs2/')

option_list = list(
  make_option(c("-c", "--cancer_type"), action="store", default="COADREAD", type='character',
              help="chose the cancer type")
)

opt = parse_args(OptionParser(option_list=option_list))
ct = opt$cancer_type   
print(ct)

# load windows and make a genomic ranges object
list_windows = list.files(path = "../RMRSigs/results/windows/window_tri_matched/800k", pattern = "*", full.names = T)
window = dplyr::bind_rows(lapply(list_windows, read_csv))
gr_window = GRanges(seqnames = window$chr, ranges = IRanges(start = window$start, end = window$end),
                    window = window$window, nt_at_risk = window$nt_at_risk)
gr_window = gr_window[seqnames(gr_window) != "chrY",] # remove chrY (because only in men)
remove(window)
sum(width(gr_window))
length(unique(gr_window$window))
gc()

# load alignability mask
alig = read_csv("../PHD_PROJECT/data_v2/filters/processed_alignability_LIFTOVERCUPs_BLACKLIST_crg75_without_exon_ctcf_ets_apobec_hairpins.csv")
gr_alig = GRanges(seqnames = alig$chr, ranges = IRanges(start = alig$start,
                                                        end = alig$end))
gr_alig_new = reduce(intersect(gr_window, gr_alig))
sum(width(gr_alig))
sum(width(gr_alig_new))
gbs = sum(width(gr_alig_new))

list_gr_window = split(gr_window, gr_window$window)
remove(gr_window)
remove(gr_alig)
remove(alig)
gc()

# read mutations 
# 1 file per sample
files = list.files('/g/strcombio/fsupek_data/users/msalvadores/RMD_signatures_project/WGS/tmp_muts', '*', full.names = T)
sample_ids = sub("\\.csv", "", files)
sample_ids = sub(".*muts_pass_([0-9a-zA-Z\\-]*)", "\\1", sample_ids)

meta = read_csv("data/metadata/comb_metadata_final_6datasets_updatedNov2021.csv",
                col_types = cols(icgc_donor_id = "c", case_submitter_id= "c"))
#meta = meta[meta$TMB >= round(3*gbs/10^6),] # 3 muts/Mb
meta = meta[grepl(ct, meta$tissue),]
sel_files = files[sample_ids %in% unique(meta$sample_id)]

cl = makeCluster(2, outfile = "")
clusterExport(cl, varlist = c("gr_alig_new", "list_gr_window", "ct", "gbs"))

RMDsplit = parLapply(cl, sel_files, function(f){
  library(GenomicRanges)
  library(readr)
  
  s = sub("\\.csv", "", f)
  s = sub(".*muts_pass_([0-9a-zA-Z\\-]*)", "\\1", s)
  print(f)
  print(s)
  
  mut <-read_csv(f)
  mut = mut[(mut$REF %in% c('A', 'T', 'C', 'G') & mut$ALT %in% c('A', 'T', 'C', 'G')),]
  gr_muts <- GRanges(seqnames = mut$chr, ranges = IRanges(start = mut$start, end = mut$end))
  #seqlevelsStyle(gr_muts) <- "UCSC"
  
  ovr = findOverlaps(gr_muts, gr_alig_new)
  gr_muts_fil = gr_muts[queryHits(ovr)]
  length(gr_muts)
  length(gr_muts_fil)
  
  tmb = length(gr_muts_fil)
  #t = sel_tmb*length(unique(gr_window$window))
  if(tmb >= round(3*gbs/10^6)){ 
    
    if(tmb > round(20*gbs/10^6)){
      gr_muts_fil = gr_muts_fil[sample(1:length(gr_muts_fil), 
                                       round(20*gbs/10^6), replace = F)]
    }
    
    count = dplyr::bind_rows(lapply(list_gr_window, function(d){
      c <- countOverlaps(d, gr_muts_fil[as.character(seqnames(gr_muts_fil)) == as.character(unique(seqnames(d)))])
      df_count = cbind(data.frame('sample_id' = s, 'cancer_type'= ct, 
                                  "count" = sum(c), "window" = unique(d$window),
                                  "nt_at_risk" = unique(d$nt_at_risk),
                                  "tmb" = tmb, 
                                  stringsAsFactors = FALSE))
      return(df_count)
    }))
    return(count)
  } # if
})# rmd split

stopCluster(cl)

# merge the split
RMD_merge = dplyr::bind_rows(RMDsplit)
dim(RMD_merge) 
print(RMD_merge[1:5,1:5])

# 4) save the file with RMD and the file with windows
dim(RMD_merge)
write.csv(RMD_merge, paste0("results/rmd_counts_updateNov2021_v2/", ct, ".csv"), row.names = FALSE)
