library(readr)
library(viridis)

wi = read_csv("data/lists/all_windows_list.csv")
head(wi)

res = dplyr::bind_rows(lapply( c(0.1, 0.2, 0.5), function(winds_affected){
  res1 = dplyr::bind_rows(lapply( c(2, 3, 5), function(increase){
    sig = paste0("sig_", winds_affected, "_", increase)
    d = wi[,c("window", "index", "chr", "start")]
    d$sig = sig
    d$value = 1
    random_wi = sample(x = 1:nrow(d),size = round(winds_affected*nrow(d)))
    d$value[random_wi] = increase
    return(d)
  }))
  return(res1)
}))

head(res)
table(res$sig)
table(res$value)
write.csv(res, "results/simulated_data/artificial_signatures_n9.csv", row.names = F)

#####
res = read_csv("results/simulated_data/artificial_signatures_n9.csv")
ggplot(res[res$chr == "chr1" & res$index < 50,], aes(x = index, y = value, col = sig)) +
  geom_line() +
  facet_wrap(vars(sig), scales = "free_x") +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(legend.position = "none")

table(res$value)
