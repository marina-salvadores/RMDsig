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
library(Biostrings)

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs/")

# 1) Get trint contexts
ms = read_csv('../PHD_PROJECT/data/genomic_features/MS96/MS_tri_6287WGS.csv', 
              col_types = cols(cancer_subtype = col_character()))
ms1 = data.frame('MS' = colnames(ms)[6:ncol(ms)], stringsAsFactors = FALSE)
remove(ms)
ms1$context = sub("(.)\\((.)>.\\)(.)", "\\1\\2\\3", ms1$MS)
ms1 = ms1[!duplicated(ms1$context),]
ms1 = bind_rows(
  lapply(split(ms1, ms1$context), function(d){
    d$context_op = as.character(reverseComplement(x = DNAStringSet(d$context)))
    return(d)
  })
)
ms1$MS = NULL

# 3) Create windows matrix
df_window_merge = read_csv("../PHD_PROJECT/data_v2/lists/all_windows/all_windows_list.csv")
head(df_window_merge)
df_window_merge$len = df_window_merge$end - df_window_merge$start +1
df_window_merge$end[df_window_merge$window == "chr19_58500002_59128984"] = df_window_merge$end[df_window_merge$window == "chr19_58500002_59128984"] -1
df_window_merge$end[df_window_merge$window == "chr4_190400002_191154277"] = df_window_merge$end[df_window_merge$window == "chr4_190400002_191154277"] -1
df_window_merge = df_window_merge[df_window_merge$len > 400000,]
min(df_window_merge$len)
table(df_window_merge$len == 10^6)

# load alignability
#### add new align here:
alig = read_csv("../PHD_PROJECT/data_v2/filters/processed_alignability_LIFTOVERCUPs_BLACKLIST_crg75_without_exon_ctcf_ets_apobec_hairpins.csv")
gr_alig = GenomicRanges::GRanges(seqnames=alig$chr, ranges = IRanges(start = alig$start, end = alig$end))

# loop here:
genes = df_window_merge$window

cl = makeCluster(40, outfile = "")
#cl = makeCluster(1, outfile = "")
clusterExport(cl, varlist = c("df_window_merge", "gr_alig", "ms1"))

per_gene = parLapply(cl, genes, function(g){
#per_gene = lapply(genes, function(g){
  
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  print(g)
  d = df_window_merge[df_window_merge$window == g,]
  gr_d = GRanges(seqnames=d$chr, ranges = IRanges(start = d$start, end = d$end))
  file_filtered = intersect(gr_d, gr_alig)
  file_filtered$len = width(file_filtered)
  nt_with_alig1 = sum(file_filtered$len)
  
  if(nt_with_alig1 >= 1000){
    # get sequence
    seq = getSeq(Hsapiens, file_filtered)
    
    # get chr frequencies
    tri = trinucleotideFrequency(x =seq)
    newtri = colSums(tri)
    df_tri = as.data.frame(newtri)
    colnames(df_tri) = "count"
    df_tri$tri = rownames(df_tri)  
    
    df_tri_fin = dplyr::bind_rows(lapply(split(ms1, ms1$context), function(d){
      a = df_tri[df_tri$tri %in% c(d$context, d$context_op),]
      new = data.frame(context = d$context, count = sum(a$count), stringsAsFactors = F)
      return(new)
    }))
    
    fin = t(df_tri_fin[,2:ncol(df_tri_fin)])
    colnames(fin) = df_tri_fin$context
    fin1 = cbind(data.frame(window = g, nt_alig1 = nt_with_alig1, stringsAsFactors = F), fin)
    return(fin1)
  }
  
}) # per gene

stopCluster(cl)

df_per_gene = dplyr::bind_rows(per_gene)
write.csv(df_per_gene, paste0("results/windows/trint_composition_3062original_windows.csv"), row.names = F)
