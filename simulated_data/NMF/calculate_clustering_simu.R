## Cluster the NMF factors, and find combination that yields best silhouette (and possibly deviation from random data)
## Also run in parallel on as many CPUs as specified
## Returns a data.table with nFact, k, and the quality measures per nFact*k combination.
#' @param maxNumClusters Global override - holds true for any number of nFact... note that some nFact will have less clusters than that 
#'                       For a given nFact, the max # clusters is 5*nFact.  Meaning: max 10 clusters for 2 factors, max 15 for 3 factors, etc.
#'                       (probably even this many does not make sense to extract...)

#initialize cosine distance function which is used for clustering
cosineDist <- function(x) {
  mat = as.dist(t(1 - cosine(t(x))))
  return(mat)
}

library(lsa)
library(readr)
library(stringr)
library(data.table)
library(dplyr)
library(cluster)
#library(clusterCrit)

#setwd("/g/strcombio/fsupek_home/msalvadores/Documents/PHD_PROJECT")
setwd("/slgpfs/projects/irb57/msalvadores/RMDSigs")

# optparse
suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
option_list = list(
  make_option(c("-c", "--cancer_type"), action="store", default="BLCA", type='character',
              help="select cancer type")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

ct = opt$cancer_type
print(ct)

# read data
nmfHmatAllByFact = readRDS(paste0("results/sig_matrix_simu_withct/", ct, "/NMF_Hmatrix.Rda"))
nmfWmatAllByFact = readRDS(paste0("results/sig_matrix_simu_withct/", ct, "/NMF_Wmatrix.Rda"))
#dataM = read_csv(paste0("results/sig_matrix_", type, "/", folder, "/df_coefs_bin_w1_", type, ".csv"))

# function 1
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y = y[!is.na(y)]
  return(y)
}

# function 2
calc_autocorr = function(mat){
  
  re = unlist(apply(mat, 1, function(d){
    newd = remove_outliers(d)
    r = acf(d, lag.max = 1, plot = FALSE)
    r = r$acf[2]
    
    r_ran = unlist(
      lapply(1:100, function(i){
        d_ran = sample(d)
        r_ran = acf(d_ran, lag.max = 1, plot = FALSE)
        return(r_ran$acf[2])
      })
    )
    
    zsco = (r - mean(r_ran)) / sd(r_ran)
    return(zsco)
  })
  )
  return(re)
} 

# function 3
#initialize function for clustering
optimizeClustersPar = function(nmfHmatAllByFact = NULL, maxNumClusters = maxK, threads = nCPUs) {
  # the highest nFact (number of factors) we had checked in this dataset... the lowest is always assumed to be ==2 but this could be inspected as well
  maxNumFact = names(nmfHmatAllByFact) %>% str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer %>% max
  
  library(parallel)
  cl = makeCluster(threads, useXDR = FALSE);
  clusterEvalQ(cl, library(cluster)) # to be able to run "pam"  (the library 'cluster' refers to data clustering, not compute-clusters!)  
  sigDiscoScores = data.table()
  
  for (aNmfSettings in names(nmfHmatAllByFact)) {
    
    # aNmfSettings is a string like "nFact=002", which encodes both the "nFact" parameter
    nFact = aNmfSettings %>% str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer
    print(nFact)
    
    # how many iterations of NMF were run for this maxNumFact.
    numNmfIters = nmfHmatAllByFact[[aNmfSettings]]$id %>% str_match("^run=(\\d+)_") %>% .[, 2] %>% as.integer %>% max
    
    nmfReBootstrappedDists = cosineDist(nmfHmatAllByFact[[aNmfSettings]][,2:ncol(nmfHmatAllByFact[[aNmfSettings]])])
    
    # parallel code: in each thread, runs the "pam" clustering for a different # clusters, starting from the same matrix of cosine similarities
    # (note: Sanger (Alexandrov et al.) always run this with #clusters = #NMF factors, which may not be ideal)
    clusterExport(cl = cl, list("nmfReBootstrappedDists"), envir = environment()) # , envir=tand)
    pamOutputByIter = parLapplyLB(# load balancing - use parLapplyLB
      cl = cl,
      X = as.list(rev(2:min(maxNumClusters, nFact * 5))), # pam() with different values of k is run in parallel 
      #X = as.list(rev(5:7)), # for proofs
      fun = function(x) {
        set.seed(42 + x);
        pam(nmfReBootstrappedDists, k = x, diss = T)
      })
    
    for (clusters in pamOutputByIter) {
      k = length(clusters$medoids);
      
      wMat <- nmfWmatAllByFact[[aNmfSettings]]
      hMat <- nmfHmatAllByFact[[aNmfSettings]]
      wMat <- t(wMat[rownames(wMat) %in% clusters$medoids, 2:ncol(wMat)])
      wMat <- as.matrix(wMat)
      hMat <- hMat[rownames(hMat) %in% clusters$medoids, 2:ncol(hMat)]
      hMat <- as.matrix(hMat)
      mult = wMat %*% hMat
      rownames(mult) = sub("^X", "", rownames(mult))
      rownames(mult) = gsub("\\.", "-", rownames(mult))
      
      sigDiscoScores = sigDiscoScores %>% rbind(data.table(nFact = nFact, k = k,
                                                           avgClScore = mean(clusters$silinfo$clus.avg.widths), minClScore = min(clusters$silinfo$clus.avg.widths),
                                                           secondWorstClScore = clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                                                           avgPtScore = mean(silhouette(clusters)[, "sil_width"]), medianPtScore = median(silhouette(clusters)[, "sil_width"])
                                                           ))
      
    }# for loop k
  } # for loop aNmfSettings
  
  stopCluster(cl)
  #rm(maxNumClusters, threads, cl, pamOutputByIter, k )
  return(sigDiscoScores)
}

#run clustering
sigClusterScores = optimizeClustersPar(nmfHmatAllByFact, maxNumClusters = 25, threads = 13)
saveRDS(sigClusterScores, paste0("results/sig_matrix_simu_withct/",ct,"/NMF_sigClusterScores.Rda")) #'quality' of clustering(s), used to determine the number of signatures 

