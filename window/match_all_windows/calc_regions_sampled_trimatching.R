suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-c", "--chr"), action="store", default="chr1", type='character',
              help="chr")
)

opt = parse_args(OptionParser(option_list=option_list))

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

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs")

# 0) Get window
#w = "chr1"
chr = opt$chr
print(chr)

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

# 2) Create windows matrix
df_window_merge = read_csv("../PHD_PROJECT/data_v2/lists/all_windows/all_windows_list.csv")
df_window_merge$end[df_window_merge$window == "chr19_58500002_59128984"] = df_window_merge$end[df_window_merge$window == "chr19_58500002_59128984"] -1
df_window_merge$end[df_window_merge$window == "chr4_190400002_191154277"] = df_window_merge$end[df_window_merge$window == "chr4_190400002_191154277"] -1

# select genes within chr 
df_window_merge = df_window_merge[df_window_merge$chr == chr,]

# load alignability
#### add new align here:
alig = read_csv("../PHD_PROJECT/data_v2/filters/processed_alignability_LIFTOVERCUPs_BLACKLIST_crg75_without_exon_ctcf_ets_apobec_hairpins.csv")
gr_alig = GenomicRanges::GRanges(seqnames=alig$chr, ranges = IRanges(start = alig$start, end = alig$end))

# 3) Get chr frequencies
tri = read_csv("results/windows/tri_composition/tri_composition_original_windows_imperfect_matching_8e+05.csv")
tri = tri[tri$window %in% df_window_merge$window,]
df_window_merge = df_window_merge[df_window_merge$window %in% tri$window,]

# get sequence
seq_chr = unmasked(Hsapiens[[chr]])

# loop here:
genes = df_window_merge$window

cl = makeCluster(20, outfile = "")
clusterExport(cl, varlist = c("df_window_merge", "gr_alig", "tri", "ms1", "chr", "seq_chr"))

per_gene = parLapply(cl, genes, function(g){
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  print(g)
  d = df_window_merge[df_window_merge$window == g,]
  gr_d = GRanges(seqnames=d$chr, ranges = IRanges(start = d$start, end = d$end))
  file_filtered = intersect(gr_d, gr_alig)
  file_filtered$len = width(file_filtered)
  nt_with_alig1 = sum(file_filtered$len)
  
  if (length(file_filtered) > 0 & nt_with_alig1 >= 500000){ #>20% of alig=1
    
    # extend matrix to 1 nt width
    ext_window = unlist(tile(file_filtered, width = 1))
    
    # calculate context
    ext_window$st = start(ext_window)-1
    ext_window$en = end(ext_window)+1
    ext_window = ext_window[ext_window$st >= d$start]
    ext_window = ext_window[ext_window$en <= d$end]
    
    seq_tmp = DNAStringSet(seq_chr, start=ext_window$st, end=ext_window$en)
    ext_window$context = as.character(seq_tmp)
    ext_window = ext_window[!grepl("N", ext_window$context)]
    
    if(length(ext_window) > 500000){
      # match context and context op
      for (con in ms1$context){
        op = ms1$context_op[ms1$context == con]
        ext_window$context[ext_window$context == op] = con
      }
      ext_window$id = paste0(seqnames(ext_window), "_", start(ext_window))
      
      # Apply match
      probs = tri[tri$window == g,]
      df_matchi = data.frame(tri = colnames(probs)[3:ncol(probs)], 
                             count = as.numeric(probs[1,3:ncol(probs)]), stringsAsFactors = F)
      
      if(!is.null(df_matchi)){
        
        res2 = lapply(split(df_matchi, df_matchi$tri), function(p){
          s = sample(x = ext_window$id[ext_window$context == as.character(p$tri)],size = p$count, replace = F)
          return(s)
        })
        
        subset = ext_window[ext_window$id %in% unlist(res2),]
        new = reduce(subset)
        
        # create matrix
        df_alig1 = data.frame("window" = d$window,
                              'chr' = seqnames(new),
                              'start' = start(new),'end' = end(new),
                              stringsAsFactors = FALSE)
        df_alig1$nt_at_risk = sum(width(new))
        return(df_alig1)
      }# if min
    }# if NNN ext_window
  }
  
}) # per gene

stopCluster(cl)

df_per_gene = dplyr::bind_rows(per_gene)
write.csv(df_per_gene, paste0("results/windows/window_tri_matched/800k/", chr, ".csv"), row.names = F)