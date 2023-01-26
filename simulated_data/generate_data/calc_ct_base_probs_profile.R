library(readr)
library(dplyr)

setwd('/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs/')

wi = read_csv("../PHD_PROJECT/data_v2/lists/all_windows/all_windows_list.csv")
meta = read_csv("../PHD_PROJECT/data_v2/metadata/metadata_complete_alig75_6datasets.csv")
meta = meta[,1:6]
table(meta$source)
table(meta$source, meta$CommonCTCode)
ovca = c(meta$sample_id[meta$source == "OVCARE"], "P00681")

#  Calculate ct-based profile 
files = list.files("results/rmd_counts_originalwindows_max", "*", full.names = T, recursive = T)

res = dplyr::bind_rows(lapply(files, function(f){
  
  print(f)
  rmd = read_csv(f)
  rmd = rmd[! rmd$sample_id %in% ovca,]

  rmd$RMD = rmd$count*10^6/rmd$nt_at_risk
  rmd1 = merge(wi[,c("window", "arm")], rmd, by = "window", all = F)
  rmd2 = dplyr::bind_rows(lapply(split(rmd1, rmd1$arm), function(d){
    a = tidyr::spread(d[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
    #a[3:ncol(a)] = a[3:ncol(a)]/rowMedians(as.matrix(a[3:ncol(a)]))
    a[3:ncol(a)] = a[3:ncol(a)]/rowMeans(as.matrix(a[3:ncol(a)]))
    b = tidyr::gather(a, window, RMD, -c(sample_id, cancer_type))
    return(b)
  }))
  
  newrmd = tidyr::spread(rmd2[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
  #rowMeans(newrmd[3:ncol(newrmd)])
  newrmd[1:5,1:5]
  
  av_rmd = colMeans(newrmd[,3:ncol(newrmd)])
  #av_rmd = av_rmd/sum(av_rmd)
  final = data.frame(window = colnames(newrmd[,3:ncol(newrmd)]),
                     av_rmd = av_rmd, stringsAsFactors = F)
  final$cancer_type = unique(newrmd$cancer_type)
  return(final)
  
}))

dim(res)
head(res)
sum(res$av_rmd[res$cancer_type == "BTCA"])
table(is.na(res$av_rmd))
write.csv(res, "results/simulated_data/ct_base_profile.csv", row.names = FALSE)



####################################
library(readr)
new = read_csv("results/simulated_data/ct_base_profile_prue.csv")
dim(new)
head(new)

head(new)
wi = read_csv("data/lists/all_windows_list.csv")

m = merge(wi, new, by = "window", all = F)
head(m)

ggplot(m[m$chr == "chr1",], aes(x = start, y = av_rmd, col = cancer_type)) +
  geom_line() +
  facet_wrap(vars(cancer_type), scales = "free_x") +
  theme_minimal() +
  theme(legend.position = "none")
