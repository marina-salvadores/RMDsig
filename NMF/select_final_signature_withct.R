library(readr)
library(dplyr)
library(cluster)
library(lsa)
suppressPackageStartupMessages(require(optparse))

setwd("/slgpfs/projects/irb57/msalvadores/RMDSigs_update_v2")

#initialize cosine distance function which is used for clustering
cosineDist <- function(x) {
  mat = as.dist(t(1 - cosine(t(x))))
  return(mat)
}

option_list = list(
  make_option(c("-c", "--cancer_type"), action="store", default="COADREAD", type='character',
              help="chose the cancer type"),
  make_option(c("-n", "--n"), action="store", default=NA, type='double',
              help="chose the cancer type"),
  make_option(c("-k", "--k"), action="store", default=NA, type='double',
              help="chose the cancer type")
)

opt = parse_args(OptionParser(option_list=option_list))
ct = opt$cancer_type   


nmfHmatAllByFact = readRDS(paste0("results/sig_matrix/NMF_Hmatrix.Rda"))
nmfWmatAllByFact = readRDS(paste0("results/sig_matrix/NMF_Wmatrix.Rda"))

#extract signatures 
# (i.e, cluster medoids) for the desired combination of 
# the number of NMF factors (nFact) and number of clusters (k)
nFact = opt$n
k = opt$k #the final number of signatures 

wMatHere = nmfWmatAllByFact[[sprintf("nFact=%03d", nFact)]]
hMatHere = nmfHmatAllByFact[[sprintf("nFact=%03d", nFact)]] # the NMF w-matrices had been transposed, so rows in the H-matrices and in the W-matrices are the same

# need to run the clustering again...
set.seed(42 + k);
# this should reproduce the same output as in the "triSig.optimizeClustersPar"
clustHere = pam(cosineDist(hMatHere[,2:ncol(hMatHere)]), k = k, diss = T)
signatures = hMatHere[clustHere$medoids,]
samps = wMatHere[clustHere$medoids,]

#plot signatures
library(ggplot2)
signatures[,1:5]

d = tidyr::gather(signatures, window, value, -c(id))
d$chr = sub("(chr[0-9XY]*)_[0-9]*_[0-9]*$", "\\1", d$window)
d$pos = as.numeric(sub("(chr[0-9XY]*)_([0-9]*)_[0-9]*$", "\\2", d$window))
head(d)
unique(d$chr)
ids = unique(d$id)

# save signatures
df_pc = as.data.frame(t(samps[,2:ncol(samps)]))
colnames(df_pc) = samps$id
df_pc$sample_id = rownames(df_pc)
df_pc$sample_id = sub("^X", "", df_pc$sample_id )
head(df_pc)

write.csv(df_pc, paste0("results/Sigs_NMF/df_pc_", k, "_", nFact, ".csv"), row.names = F)

head(d)
write.csv(d, paste0("results/Sigs_NMF/signatures_", k, "_", nFact, ".csv"), row.names = F)

