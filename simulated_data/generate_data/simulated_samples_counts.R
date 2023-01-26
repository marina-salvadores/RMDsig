library(readr)
library(dplyr)
library(sampling)
library(parallel)
setwd("/g/strcombio/fsupek_home/msalvadores/Documents/RMRSigs/")

# load samples
samples = read_csv("data/metadata/metadata_samples_rmd_max.csv")
meta = read_csv("../PHD_PROJECT/data_v2/metadata/metadata_complete_alig75_6datasets.csv")
ovca = c(meta$sample_id[meta$source == "OVCARE"], "P00681")
samples = samples[!samples$sample_id %in% ovca,]
head(samples)

wi = read_csv("data/lists/all_windows_list.csv")
head(wi)

ct_base = read_csv("results/simulated_data/ct_base_profile.csv")
ct_base = ct_base[order(ct_base$window),]
sigs = read_csv("results/simulated_data/artificial_signatures_n9.csv")
sigs = sigs[order(sigs$window),]

head(ct_base)
head(sigs)

for(ct in unique(ct_base$cancer_type)){

  print(ct)
  df_samples = samples[samples$cancer_type == ct,]
  ct_ind = ct_base[ct_base$cancer_type == ct, ]
  
  for(sig_contrib in c(0.1, 0.2, 0.4)){
    print(sig_contrib)
    
    for(samp_affected in c(0.05, 0.1, 0.2)){
      print(samp_affected)
      
      df_samples_affected_distrib = dplyr::bind_rows(
        lapply(unique(sigs$sig), function(s){
          new = df_samples[,c("sample_id", "cancer_type", "final_tmb")]
          new$sig = s
          new$sample_affected = 0
          new$sample_affected[sample(1:nrow(new), size = round(nrow(new)*samp_affected))] = 1
          return(new)
        })
      )
      freq = data.frame(table(df_samples_affected_distrib$sample_id, df_samples_affected_distrib$sample_affected))
      freq = freq[freq$Var2 == 1,]
      table(freq$Freq > 2)
      
      if ( any(freq$Freq > 2) ) {
        sel = freq[freq$Freq > 2,]
        s1 = as.character(sel$Var1[1])
        for(s1 in sel$Var1){
          num = sel$Freq[sel$Var1 == s1]
          df_samples_affected_distrib$sample_affected[df_samples_affected_distrib$sample_id == s1 & 
          df_samples_affected_distrib$sample_affected == 1]  = sample(c(rep(1,2), rep(0, num-2)), num)
        }
      } # if
      
      df_samples_affected_distrib$samp_affected = samp_affected
      df_samples_affected_distrib$sig_contribution = sig_contrib
      name = paste0("/g/strcombio/fsupek_data/users/msalvadores/RMD_signatures_project/simulated_data/index/index_",
                    ct, "_sigcontrib_", sig_contrib, "_samps_", samp_affected, ".csv")
      write.csv(df_samples_affected_distrib, name, row.names = F)
      
      cl = makeCluster(40, outfile = "")
      clusterExport(cl, varlist = c("ct_ind", "sigs"))

      fin = parLapply(cl = cl, X = split(df_samples_affected_distrib,
                      df_samples_affected_distrib$sample_id), fun = function(d){
        
        library(sampling)
        tmb = unique(d$final_tmb)
        ct = unique(d$cancer_type)
        sig_contrib = unique(d$sig_contribution)
        d$per = d$sample_affected*sig_contrib
        per_ct = 1 - sum(d$per)
        red_d = d[d$sample_affected == 1,]
        
        if(nrow(red_d) > 0){

          pik_ct =ct_ind[,c("window", "av_rmd")]
          pik_ct$sig = "ct"
          pik_ct = tidyr::spread(pik_ct[,c("window", "sig", "av_rmd")], window, av_rmd)
          pik_ct[,2:ncol(pik_ct)] = pik_ct[,2:ncol(pik_ct)]/rowSums(pik_ct[,2:ncol(pik_ct)])
          pik_ct[,2:ncol(pik_ct)] = per_ct*tmb*pik_ct[,2:ncol(pik_ct)]
          
          n = dplyr::bind_rows(lapply(unique(red_d$sig), function(s){
            per_sig = red_d$per[red_d$sig == s]
            prof = merge(x = ct_ind, y = sigs[sigs$sig == s,], 
                         by = "window", all = F)
            prof$prof = prof$av_rmd*prof$value
            df_prof = tidyr::spread(prof[,c("window", "sig", "prof")], window, prof)
            df_prof[,2:ncol(df_prof)] = df_prof[,2:ncol(df_prof)]/rowSums(df_prof[,2:ncol(df_prof)])
            df_prof[,2:ncol(df_prof)] = per_sig*tmb*df_prof[,2:ncol(df_prof)]
            return(df_prof)
          }))
          
          table(colnames(n) == colnames(pik_ct))
          n2 = dplyr::bind_rows(n, pik_ct)
          final_pik = colSums(n2[,2:ncol(n2)])
          simu = UPmultinomial(pik = final_pik)
          df_simu = as.data.frame(simu)
          df_simu$window = colnames(n2[,2:ncol(n2)])

        }else{
          pik_ct =ct_ind[,c("window", "av_rmd")]
          pik_ct$sig = "ct"
          pik_ct = tidyr::spread(pik_ct[,c("window", "sig", "av_rmd")], window, av_rmd)
          pik_ct[,2:ncol(pik_ct)] = pik_ct[,2:ncol(pik_ct)]/rowSums(pik_ct[,2:ncol(pik_ct)])
          pik_ct[,2:ncol(pik_ct)] = per_ct*tmb*pik_ct[,2:ncol(pik_ct)]
          
          final_pik = pik_ct[,2:ncol(pik_ct)]
          simu = UPmultinomial(pik = final_pik)
          df_simu = as.data.frame(simu)
          df_simu$window = colnames(pik_ct[,2:ncol(pik_ct)])
        } # if
        
        head(df_simu)
        df_simu$tmb = tmb
        df_simu$sample_id = unique(d$sample_id)
        df_simu$cancer_type = unique(d$cancer_type)

        return(df_simu)
      })
      df_fin = dplyr::bind_rows(fin)
      stopCluster(cl)
      
      name_fin = paste0("/g/strcombio/fsupek_data/users/msalvadores/RMD_signatures_project/simulated_data/simu_counts/simu_counts_",
                    ct, "_sigcontrib_", sig_contrib, "_samps_", samp_affected, ".csv")
      write.csv(df_fin, name_fin, row.names = F)
      
    } 
  } 
}

