# this script divides the genome into 1Mb exact windows, filtering out all regions
# outside alignability75 mask + other removed regions

library(readr)
library(GenomicRanges)
library(parallel)

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs")

# select window size
bin = 1*10^6

#get chr arm coordinates
print(bin)
cy = read_tsv('../PHD_PROJECT/data/supp/table_bowser_cytoBand') # 0 based
cy$window = paste0(cy$`#chrom`, "_", cy$name)
cy$arm = sub("(chr[0-9X]*\\_[pq]+).*","\\1",cy$window)
gr_cy = GRanges(seqnames = cy$`#chrom`, ranges = IRanges(start = cy$chromStart+1, 
                end = cy$chromEnd), name = cy$name, window = cy$window, arm = cy$arm)
gr_cy$len = width(gr_cy)

# alignability 75 + remove problematic parts (exons, ctcf peaks, etc...)
alig = read_csv("../PHD_PROJECT/data_v2/filters/processed_alignability_LIFTOVERCUPs_BLACKLIST_crg75_without_exon_ctcf_ets_apobec_hairpins.csv")
gr_alig = GRanges(seqnames = alig$chr, ranges = IRanges(start = alig$start, end = alig$end))
sum(width(gr_alig))

#by cytoband
cytos = unique(cy$arm)

# paralelize 
cl = makeCluster(20, outfile = "")
clusterExport(cl, varlist = c("gr_cy", "gr_alig", "bin"))

final = parLapply(cl, cytos, function(cyt){
  library(GenomicRanges)
  
  print(cyt)
  gr = gr_cy[gr_cy$arm == cyt] # get coordinates for the chr arm
  gr_inter = intersect(gr, gr_alig) # keep only those that intersect with the alignability mask
  gr_inter = reduce(gr_inter)
  nts = sum(width(gr_inter)) # total number of usable nts in the chr arm
  run = trunc(nts/bin) # times 1 mb is in total nts
  max = trunc(nts/run) # when divide total nts into the num of runs
  print(run)
  print(max)
  
  if(run > 0){
    ext_window = unlist(tile(gr_inter, width = 1))
    
    res = lapply(1:run, function(j){ # run 1 per window
      print(j)
      a = (j-1)*max # coordinate start window
      b = (j)*(max-1) # coordinate end window
      
      gr_sel = ext_window[a:b]
      gr_sel = gr_sel[sample(1:length(gr_sel), bin)] # from the window
      # sample 1Mb positions to have 1 mb exact window
      gr_window = reduce(gr_sel)
      print(sum(width(gr_window))) # width = 1Mb
      window = data.frame(window = paste0(cyt, "_w", j), "ntrisk1" = sum(width(gr_window)),
                          'chr' = as.character(seqnames(gr_window)),
                          'start' = start(gr_window),'end' = end(gr_window),
                          stringsAsFactors = FALSE)
      return(window)
    })
    
    window_fin = dplyr::bind_rows(res)
    return(window_fin)
  }#if
  
})

stopCluster(cl)

df_fin = dplyr::bind_rows(final)
print(length(unique(df_fin$window)))
print(dim(df_fin))
print(df_fin[1:5,])
# save the results
write.csv(df_fin, "results/windows/windows_1mb_per_arm.csv", row.names = F)
