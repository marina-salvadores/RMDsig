library(readr)

a = read_csv("results/windows/tri_composition/tri_composition_original_windows_imperfect_matching_8e+05.csv")
a[1:5,1:5]

a = a[,3:ncol(a)]
b = a/rowSums(a)
rowSums(b)
b[1:5,1:5]
write.csv(b, "results/windows/tri_composition/tri_composition_original_windows_imperfect_matching_8e+05_feqmatrix.csv", row.names = F)


head(b)
library(matrixStats)
hist(colSds(as.matrix(b)))

range(b)
hist(colRanges(as.matrix(b)), breaks = 100)
