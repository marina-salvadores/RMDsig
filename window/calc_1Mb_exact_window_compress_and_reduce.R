library(readr)

setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs")

df_fin = read_csv("results/windows/windows_1mb_per_arm.csv")
head(df_fin)
length(unique(df_fin$window))

res = lapply(split(df_fin, df_fin$window), function(d){
  new = d[1,]
  new$start = min(d$start)
  new$end = max(d$end)
  return(new)
})

df_res = dplyr::bind_rows(res)

df_res$len = df_res$end - df_res$start +1
hist(df_res$len, breaks = 1000)
hist(df_res$len, breaks = 1000, xlim = c(0, 2000000))
table(df_res$len < 3*10^6 & df_res$len > 1*10^6)

df_res1 = df_res[df_res$len < 3*10^6,]
hist(df_res$len)
hist(df_res1$len, breaks = 100)

length(unique(df_res1$window))
write.csv(df_res1, "results/windows/windows_1mb_per_arm_compress.csv", row.names = F)

newdf_fin = df_fin[df_fin$window %in% df_res1$window, ]
length(unique(df_fin$window))
length(unique(newdf_fin$window))
write.csv(newdf_fin, "results/windows/windows_1mb_per_arm_2390reduce.csv", row.names = F)


#######################################################
# profiles
me = rmd
head(me)
win = read_csv("results2/windows_1mb_per_arm_compress.csv")
me = merge(me, win, by ="window", all = F)
head(me)
table(me$arm)
#me = me[me$chr != "chrY",]

ggplot(me, aes(x = start, y = count)) +
  geom_line() +
  #geom_vline(data = sub, aes(xintercept = x), linetype = "dashed") +
  facet_wrap(vars(arm), scales = "free_x") +
  #scale_color_manual(values = c("blue", "navy")) +
  theme_minimal() +
  theme(legend.position = "none")


length(unique(df_res1$window))
write.csv(df_res1, "results2/windows_1mb_per_arm_compress.csv", row.names = F)


#########################
