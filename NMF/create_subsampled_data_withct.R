library(readr)
library(dplyr)
library(optparse)
library(sampling)
library(parallel)

setwd("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2")

# optparse
suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
option_list = list(
  make_option(c("-c", "--cancer_type"), action="store", default="headneck", type='character',
              help="select cancer type")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

ct = opt$cancer_type
print(ct)


wi = read_csv("data/all_windows_list.csv")
files = list.files(path = "results/rmd_counts_updateNov2021_v2", pattern = ct, 
                   full.names = T)

# read file
rmd = dplyr::bind_rows(lapply(files, read_csv))
#table(rmd$count==0) 
table(is.na(rmd))
rmd[1:5,1:5]

tmb = rmd %>% group_by(sample_id) %>% summarise(tmb = sum(count))
#head(tmb)
#hist(tmb$tmb)
#table(tmb$tmb > 3000)
#3*1958707652*10^-6
#rmd = rmd[rmd$sample_id %in% tmb$sample_id[tmb$tmb >= 3917],]

nt = rmd[!duplicated(rmd$window),]
av_nt = trunc(mean(nt$nt_at_risk))

rmd$RMD = rmd$count*av_nt/rmd$nt_at_risk
rmd1 = merge(wi[,c("window", "arm")], rmd, by = "window", all = F)
rmd2 = dplyr::bind_rows(lapply(split(rmd1, rmd1$arm), function(d){
  a = tidyr::spread(d[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
  row_na = rowSums(is.na(a))
  a = a[row_na == 0,]
  a[3:ncol(a)] = a[3:ncol(a)]/rowMeans(as.matrix(a[3:ncol(a)]), na.rm = T)
  b = tidyr::gather(a, window, RMD, -c(sample_id, cancer_type))
  return(b)
}))

sprmd1 = tidyr::spread(rmd1[,c("sample_id", "cancer_type", "window", "RMD")], 
                       window, RMD)
sprmd1[1:5,1:5]
table(is.na(sprmd1))
row_na = rowSums(is.na(sprmd1[,2:ncol(sprmd1)]))
table(row_na == 0)
sprmd1 = sprmd1[row_na == 0,]

newrmd = tidyr::spread(rmd2[,c("sample_id", "cancer_type", "window", "RMD")], window, RMD)
newrmd[1:5,1:5]
print("NAs")
print(table(is.na(newrmd)))
row_na = rowSums(is.na(newrmd[,2:ncol(newrmd)]))
table(row_na == 0)
newrmd = newrmd[row_na == 0,]

# multiply by av RMD
newrmd[,3:ncol(newrmd)] = newrmd[,3:ncol(newrmd)]*rowMeans(sprmd1[,3:ncol(sprmd1)])
list_rmd = split(newrmd, newrmd$sample_id)

cl = makeCluster(5, outfile = "")
clusterExport(cl, varlist = c("list_rmd", "ct"))

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
  folder = paste0("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2/results/input_coefs/", ct)
  if (!file.exists(folder)){
    dir.create(file.path(folder), recursive = TRUE)
  } 
  
  write.csv(final, paste0(folder,"/coef_", inter, ".csv"), row.names = F)
})

stopCluster(cl)
