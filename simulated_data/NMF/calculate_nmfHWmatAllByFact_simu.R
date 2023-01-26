library(readr)
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


#organize and save data
# mat H
files = list.files(path = paste0("results/sig_matrix_simu_withct/", ct,"/matH"), 
                   pattern = "*", full.names = T)
nmfHmatAllByFact = list()

for(f in files){
  df = read_csv(f)
  if(nrow(df)>1){
    for(i in c(5,10,15,20,25)){
      if(i >=10){
        c = paste0("nFact=0", i)
      }else{
        c = paste0("nFact=00", i)
      }
      red = df[grepl(c, df$id),]
      nmfHmatAllByFact[[c]] = rbind(nmfHmatAllByFact[[c]], red)
    }
  }
}
names(nmfHmatAllByFact)
saveRDS(nmfHmatAllByFact, paste0("results/sig_matrix_simu_withct/", ct, "/NMF_Hmatrix.Rda"))

#mat W
files = list.files(path = paste0("results/sig_matrix_simu_withct/", ct, "/matW/"), 
                   pattern = "*", full.names = T)
nmfWmatAllByFact = list()

for(f in files){
  df = read_csv(f)
  if(nrow(df)>1){
    for(i in c(5,10,15,20,25)){
      if(i >=10){
        c = paste0("nFact=0", i)
      }else{
        c = paste0("nFact=00", i)
      }
      red = df[grepl(c, df$id),]
      nmfWmatAllByFact[[c]] = rbind(nmfWmatAllByFact[[c]], red)
    }
  }
}

names(nmfWmatAllByFact)
saveRDS(nmfWmatAllByFact, paste0("results/sig_matrix_simu_withct/", ct, "/NMF_Wmatrix.Rda")) 

