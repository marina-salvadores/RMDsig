library(readr)

files = list.files("/g/strcombio/fsupek_data/users/msalvadores/RMD_signatures_project/simulated_data/simu_counts/",
                   pattern = "*", full.names = T)

types = unique(sub(".*simu_counts_[A-Z]*\\_(sigcontrib.*)\\.csv", "\\1", files))

for (t in types) {
  print(t)
  red_files = files[grepl(t, files)]
  res = dplyr::bind_rows(lapply(red_files, function(f){
    df = read_csv(f)
    return(df)
  }))
  head(res)
  dim(res)  
  table(unique(res$cancer_type))
  write.csv(res, paste0("/g/strcombio/fsupek_data/users/msalvadores/RMD_signatures_project/simulated_data/simu_counts_pooled/",
                        t, ".csv"), row.names = F)
}
