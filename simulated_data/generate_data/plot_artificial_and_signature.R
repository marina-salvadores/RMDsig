library(readr)
library(ggplot2)
library(lsa)

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs")

files = list.files(path = "results/simulated_data/fin_simu", pattern = "signatures",full.names = T)
f = files[9]
f

sig = read_csv(f)
sig = sig[,c("id", "window", "value")]
head(sig)

real_sigs = read_csv("results/simulated_data/artificial_signatures_n9.csv")
real_sigs = real_sigs[,c("window", "sig", "value")]
head(real_sigs)

new_sigs = unique(sig$id)
real = unique(real_sigs$sig)

res = lapply(new_sigs, function(n){
  res2 = lapply(real, function(r){
    
    df_n = sig[sig$id == n,]
    df_r = real_sigs[real_sigs$sig == r,]
    
    me = merge(df_n, df_r, by = "window", all = F)
    return(data.frame(n = n, r = r, cor = cor(me$value.x, me$value.y),
                      cos = cosine(me$value.x, me$value.y), stringsAsFactors = F))
  })
  return(dplyr::bind_rows(res2))
})

res = dplyr::bind_rows(res)

head(res)
str(res)
table(res$cos > 0.80)

mat = tidyr::spread(res[,c("n", "r", "cos")], r, cos)
#mat = tidyr::spread(res[,c("n", "r", "cor")], r, cor)
mat[1:5,1:5]
rownames(mat) = mat$n
mat = mat[,2:ncol(mat)]
library(ComplexHeatmap)
Heatmap(mat,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

head(res)
red = res[!grepl("91", res$n),]
red = red[!grepl("26", red$n),]
red = red[red$cos > 0.8,]

head(real_sigs)
head(sig)

comb = merge(sig, real_sigs, by = "window", all = F)
head(comb)

comb$comb = paste0(comb$sig.x, " & ", comb$sig.y)
red$comb = paste0(red$n, " & ", red$r)

wi = read_csv("data/lists/all_windows_list.csv")
comb = merge(comb, wi, by = "window", all = F)

comb = comb[comb$comb %in% red$comb,]
head(comb)

comb = comb[comb$chr == "chr1" & comb$start < 1.5*10^8,]
ggplot() +
  geom_line(data = comb, aes(x =start, y = value.x), col = "blue") +
  geom_line(data = comb, aes(x =start, y = value.y*0.005), col = "red") +
  facet_wrap(vars(comb), scales = "free") +
  theme_minimal()
