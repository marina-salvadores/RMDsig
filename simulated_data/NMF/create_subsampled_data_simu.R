library(readr)
library(dplyr)
library(optparse)
library(sampling)
library(parallel)

setwd("/slgpfs/projects/irb57/msalvadores/RMDSigs")

option_list = list(
  make_option(c("-f", "--file"), action="store", default="NA", type='character',
              help="select file")
)

opt = parse_args(OptionParser(option_list=option_list))
file = opt$file
name = sub("\\.csv", "", file)
print(file)
print(name)

path = "results/simu_counts_pooled/"

# data required
meta = read_csv("data/metadata_complete_alig75_6datasets.csv")
#meta = read_csv("../PHD_PROJECT/data_v2/metadata/metadata_complete_alig75_6datasets.csv")
ovca = c(meta$sample_id[meta$source == "OVCARE"], "P00681") # plus 1 problematic sample from SARC
wi = read_csv("data/all_windows_list.csv")
#wi = read_csv("data/lists/all_windows_list.csv")

#read file
rmd = read_csv(paste0(path, file, ".csv"))
rmd = rmd[! rmd$sample_id %in% ovca,] # remove ovcare
rmd = rmd[rmd$window != "chr10_1_1000000",]

table(rmd$simu==0) 
rmd[1:5,1:5]
nt = read_csv("data/metadata/nt_at_risk.csv")
rmd = merge(rmd, nt, by = "window", all = F)
rmd[1:5,]

av_nt = trunc(mean(nt$nt_at_risk))
  
rmd$RMD = rmd$simu*av_nt/rmd$nt_at_risk
rmd1 = merge(wi[,c("window", "arm")], rmd, by = "window", all = F)
rmd2 = dplyr::bind_rows(lapply(split(rmd1, rmd1$arm), function(d){
  a = tidyr::spread(d[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
  a[3:ncol(a)] = a[3:ncol(a)]/rowMeans(as.matrix(a[3:ncol(a)]))
  b = tidyr::gather(a, window, RMD, -c(sample_id, cancer_type))
  return(b)
}))
  
sprmd1 = tidyr::spread(rmd1[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
sprmd1[1:5,1:5]
  
newrmd = tidyr::spread(rmd2[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
newrmd[1:5,1:5]
  
# multiply by av RMD
newrmd[,3:ncol(newrmd)] = newrmd[,3:ncol(newrmd)]*rowMeans(sprmd1[,3:ncol(sprmd1)])

list_rmd = split(newrmd, newrmd$sample_id)
  
cl = makeCluster(10, outfile = "")
clusterExport(cl, varlist = c("list_rmd", "name"))
  
parLapply(cl, 1:100, function(inter){
    
    print(inter)
    
    res1 = lapply(list_rmd, function(d){
      print(unique(d$sample_id))
      x = d[,3:ncol(d)]
      x[x<0] = 0
      d[,3:ncol(d)] = sampling::UPmultinomial(pik = x) 
      #this changes the row a bit while keeping the total sum of mutations the same
      return(d)
    })
    
    final = dplyr::bind_rows(res1)  
    print(dim(final))
    print(final[1:5,1:6])
    
    # save
    folder = paste0("/slgpfs/projects/irb57/msalvadores/RMDSigs/results/input_coefs_simu_withct/", name)
    if (!file.exists(folder)){
      dir.create(file.path(folder), recursive = TRUE)
    } 
    
    write.csv(final, paste0(folder,"/coef_", inter, ".csv"), row.names = F)
  })
  
  stopCluster(cl)
}

