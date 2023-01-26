library(readr)
library(dplyr)
library(parallel)
library(NMF)

setwd("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2")

# optparse
suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
option_list = list(
  make_option(c("-d", "--dataframe"), action="store", type = "character",
              help="file with coefs"),
  make_option(c("-f", "--minimum_factors"), action="store", default = 1,
              help="Minimum number of factors to optimiza the NMF"),
  make_option(c("-F", "--maximum_factors"), action="store", default = 20,
              help="Maximum number of factors to optimiza the NMF"),
  make_option(c("-p", "--no_cores"), action="store", default = 10,
              help="number of cores to use")
)


opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# params
file_name = opt$dataframe
ct = "pancan"
minimum_factor_range = opt$minimum_factors
maximum_factor_range = opt$maximum_factors
no_cores = opt$no_cores
factor_range = rev(seq(from = minimum_factor_range, to = maximum_factor_range, by = 1))

# Load data:
path = list.files(path = "results/input_coefs", pattern = file_name,
                  full.names = T, recursive = T)
print(path)
nIter = as.numeric(sub(".*coef_([0-9]*)\\.csv", "\\1", file_name))
print(nIter)

resamp = dplyr::bind_rows(lapply(path, read_csv))
dim(resamp)
resamp[1:5,1:5]
unique(resamp$cancer_type)
table(resamp$cancer_type)

# prepare the cluster for parallelize
cl = makeCluster(no_cores, outfile = "")
clusterExport(cl=cl, list("resamp", "nIter"),envir = environment())

total = parLapply(cl, factor_range, function(nFact){
  
  library(NMF)
  
  idString = sprintf( "nFact=%03d", nFact)
  print(idString)
  
  mat = as.matrix(resamp[,3:ncol(resamp)])
  #print(colnames(mat)[colSums(mat)==0])
  mat = mat[,colSums(mat)!=0]
  res = NMF::nmf(x = mat,rank = nFact, maxIter=10000)
  
  #recover H and W
  matH = res@fit@H
  matW = res@fit@W
  
  # change the names to the signatures/factors
  rownames(matH) = sprintf("run=%03d_%s_fact%02d", nIter, idString, 1:nFact)
  colnames(matH) <- colnames(mat)
  
  rownames(matW) = resamp$sample_id
  colnames(matW) = sprintf("run=%03d_%s_fact%02d", nIter, idString, 1:nFact)
  
  out = list()
  out[[paste0("matH_", idString)]] = matH
  out[[paste0("matW_", idString)]] = t(matW)
  return(out) 
})

stopCluster(cl)

totalH = dplyr::bind_rows(
  lapply(total, function(pair){
    sel = names(pair)[grepl("matH", names(pair))]
    a = data.frame( id = rownames(pair[[sel]]),pair[[sel]], stringsAsFactors = F)
    return(a)
  })
)

# save
folder = paste0("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2/results/sig_matrix/matH")
if (!file.exists(folder)){
  dir.create(file.path(folder), recursive = TRUE)
} 
print(totalH[1:5,1:5])
write.csv(totalH, paste0(folder, "/matH_run_", nIter, ".csv"), row.names = F)


totalW = dplyr::bind_rows(
  lapply(total, function(pair){
    sel = names(pair)[grepl("matW", names(pair))]
    a = data.frame( id = rownames(pair[[sel]]),pair[[sel]], stringsAsFactors = F)
    return(a)
  })
)

# save
folder = paste0("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2/results/sig_matrix/matW")
if (!file.exists(folder)){
  dir.create(file.path(folder), recursive = TRUE)
} 
print(totalW[1:5,1:5])
write.csv(totalW, paste0(folder, "/matW_run_", nIter, ".csv"), row.names = F)


