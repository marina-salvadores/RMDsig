# ----- helper functions
setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs/")
library(readr)

generateJitteryRow = function(freqs, numNt, jitter) {
  trunc( freqs * runif(length(freqs), min=1-jitter, max=1+jitter) * numNt)
}
# normalize each row to sum to 1. There may be a faster way.
rowNorm =function(m) {
  t( apply(m, 1, function(x) { x/sum(x) } ) )
}
euclidean = function(a, b) {
  sqrt(sum((a - b)^2))
} 

files = list.files(path = "results/windows/penta_composition", pattern = "imperfect", full.names = T)
#files  = files[!grepl("feqmatrix", files)]
#names(files) = c(100, 1000, 200, 300, 400, 500, 600, 700, 800, 900)
files  = files[!grepl("freqmatrix", files)]
names(files) = c(10, 50, 100, 1000, 200, 300, 400, 500, 600, 700, 800, 900)
files
#n = names(files)[1]
#print(n)

res = dplyr::bind_rows(lapply(names(files), function(n){
  f = files[n]
  all_counts = read_csv(f)
  all_counts = all_counts[all_counts$nt_alig1 >= 500*10^3,]
  counts = as.matrix(all_counts[,3:ncol(all_counts)])
  
  freqs=rowNorm(counts);
  meanFreqs=colMeans(freqs, na.rm = T);
  offender = which.max( apply(freqs, 1, function(x){ euclidean(x,meanFreqs) })  );
  diffs = meanFreqs - freqs[offender,];
  correctableCol = which.min(diffs)
  worstCol = which.max(abs(diffs))  # this is sometimes the same as the correctable col
  
  # note that tolerance is expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal
  tolerance = abs(diffs[worstCol]) 
  tolerance
  
  return(data.frame(iter = n, tol = tolerance, nts = nrow(counts), stringsAsFactors = F))
}))

res

library(ggplot2)
res$iter = as.numeric(res$iter)
a = ggplot(res, aes(x = iter*1000, y = tol)) +
  geom_point() +
  geom_line() +
  ylab("tolerance") +
  xlab("iterations") +
  theme_minimal()

head(res)
b = ggplot(res, aes(x = iter*1000, y = nts)) +
  geom_point() +
  geom_line() +
  ylab("windows >500k") +
  xlab("iterations") +
  theme_minimal()

library(cowplot)
plot_grid(a,b, nrow = 1)
