#Light:Dark Datasets Method Application 

library(tidyverse)
library(MetaCycle)
library(matrixStats)
library(limma)

#GRO-seq data -----
    #QC
    load("gro.rda")
    gro1 <- gro[rowSums(gro == 0) <= 4, ]
    #MC
    meta2d(infile="GRO.csv",filestyle="csv",outdir="LD_Circadian",
           minper=20, maxper=28,
           timepoints=seq(1, 22, by=3),
           outIntegration="onlyIntegration", ARSdefaultPer=24,
           outRawData=FALSE)
    #BC 
    write_tsv(gro1, "GRO.tsv")
      #Rscript BioCycle.R -i /Users/mwen/Desktop/LD_circadian/GRO.tsv -o /Users/mwen/Desktop/ -s 20  -e 28 -o /Users/mwen/Desktop/ -s 20  -e 28
    #RAIN
    library(rain)
    gro_rain <- rain(t(gro1[,-1]),deltat=3, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                     nr.series = 1,  adjp.method = "ABH", verbose = FALSE)
    #eJTK 
    write.table(gro1, "gro.txt", row.names = F, col.names = T)
    # ./eJTK-CalcP.py -f example/gro.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT
    
#Nascent-seq data ----
    #QC 
    load("Nascent.rda")
    Nascent1 <- Nascent[rowSums(Nascent == 0) <= 6, ]
    #MC 
    write.csv(Nascent1, file="Nascent1.csv", row.names=F)
    meta2d(infile="Nascent1.csv",filestyle="csv",outdir="LD_Circadian",
           minper=20, maxper=28,
           timepoints=rep(seq(0, 20, by=4),times=2),
           outIntegration="onlyIntegration", ARSdefaultPer=24,
           outRawData=FALSE)
    #BC
    library(readr)
    write_tsv(Nascent1, "Nascent.tsv")
      #Rscript BioCycle.R -i /Users/mwen/Desktop/Nascent.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
    
    #RAIN
    nascent_rain <- rain(t(Nascent1[,-1]),deltat=4, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                         nr.series = 2,  adjp.method = "ABH", verbose = FALSE)
    nascent_rain$rain_q <- qvalue(p=Nascent_rain$pVal)$qvalues 
    #eJTK
    write.table(Nascent1, "Nascent.txt", row.names = F, col.names = T)
   
#RNA-seq data -----
    load("RNA.rda")
    RNA1 <- RNA[rowSums(RNA == 0) <= 6, ]
    #Metacycle 
    library(MetaCycle)
    write.csv(RNA1, file="RNA1.csv", row.names=F)
    meta2d(infile="RNA1.csv",filestyle="csv",outdir="LD_Circadian",
           minper=20, maxper=28,
           timepoints=rep(seq(2, 22, by=4),times=2),
           outIntegration="onlyIntegration", ARSdefaultPer=24,
           outRawData=FALSE)
    #BioCycle 
    write_tsv(RNA1, "RNA.tsv")
      #Rscript BioCycle.R -i /Users/mwen/Desktop/RNA.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
    
    
    #RAIN
    library(rain)
    library("magrittr")
    for (i in 2:13) {
      RNA[,i] <- as.numeric(as.character(RNA[,i]))
    }
    
    RNA_rain <- rain(t(RNA1[,-1]),deltat=4, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                     nr.series = 2,  adjp.method = "ABH", verbose = FALSE)
    
    RNA_rain$rain_q<- qvalue(p=RNA_rain$pVal)$qvalues 
    
    #eJTK
    write.table(RNA1, "RNA.txt", row.names = F, col.names = T)
    
#XR-seq data -----
    load("XR.rda")
      #XR1 <- XR[rowSums(XR == 0) <= 6, ]
    #MetaCycle
    library(MetaCycle)
    write.csv(XR, file="XR.csv", row.names=F)
    meta2d(infile="XR.csv",filestyle="csv",outdir="LD_Circadian",
           minper=20, maxper=28,
           timepoints=rep(seq(0, 20, by=4),times=2),
           outIntegration="onlyIntegration", ARSdefaultPer=24,
           outRawData=FALSE)
    
    #BIOCycle 
    write_tsv(XR, "XR.tsv")
      #Rscript BioCycle.R -i /Users/mwen/Desktop/XR.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
    
    #RAIN
    library(rain)
    for (i in 2:13) {
      XR[,i] <- as.numeric(as.character(XR[,i]))
    }
    
    XR_rain <- rain(t(XR[,-1]),deltat=4, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                    nr.series = 2,  adjp.method = "ABH", verbose = FALSE)
    
    #eJTK
    write.table(XR, "XR.txt", row.names = F, col.names = T)
    #import results
    gro_meta2d<- read.csv(file="meta2d_GRO.csv", header=TRUE, sep=",")
    nascent_meta2d<- read.csv(file="meta2d_nascent1.csv", header=TRUE, sep=",")
    RNA_meta2d<- read.csv(file="meta2d_RNA.csv", header=TRUE, sep=",")
    XR_meta2d<- read.csv(file="meta2d_XR.csv", header=TRUE, sep=",")
    
    gro_bio <- read_tsv(file="gro_BioCycle.tsv", col_names=T) 
    nascent_bio <- read_tsv(file="Nascent_BioCycle.tsv", col_names=T) 
    RNA_bio <- read_tsv(file="RNA_BioCycle.tsv", col_names=T) 
    XR_bio <- read_tsv(file="XR_BioCycle.tsv", col_names=T) 
    
    gro_eJTK <- read.table(file="gro_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T, sep = "\t")
    nascent_eJTK <- read.table(file="Nascent_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T, sep = "\t")
    RNA_eJTK <- read.table(file="RNA_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T, sep = "\t")
    XR_eJTK <- read.table(file="XR_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T, sep = "\t")
    
    
# combine all the results -----
    combine_dataset=function(name,...){
      name1=paste(name,"bio",sep = "_")
      name2=paste(name,"rain",sep = "_")
      name3=paste(name,"meta2d",sep = "_")
      name4=paste(name, "eJTK", sep="_")
      full = cbind(cbind(cbind(get(name1)[,c(1, 2, 3)], 
                               get(name2)[,c(1,5)]),
                         get(name3)[,c(2,3,7,8,12,13)]), 
                   get(name4)[,c(19,20)])
      names(full) = c("ID", "Bio_P", "BIO_Q", "rain_P", "rain_Q", "JTK_pvalue","JTK_BH.Q",
                      "LS_pvalue","LS_BH.Q","meta2d_pvalue","meta2d_BH.Q","eJTK_pval","eJTK_BH")
      
      return(full)
    }
    
    #full datasets
    nascent_full <- combine_dataset("nascent")
    XR_full <- combine_dataset("XR")
    RNA_full <- combine_dataset("RNA")
    gro_full = gro_bio %>% 
      cbind(gro_rain) %>%
      left_join(gro_meta2d, by="ID") %>%
      left_join(gro_eJTK, by="ID") %>%
      subset(select=c(1, 2, 3,17, 21, 27, 28, 32, 33, 37, 38, 61, 62)) 
    
    save(nascent_full, file="nascent_full.rda")
    save(XR_full, file="XR_full.rda")
    save(RNA_full, file="RNA_full.rda")
    save(gro_full, file="gro_full.rda")
    
    
    

