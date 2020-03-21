library(tidyverse)

setwd("/Users/mwen/Desktop/Research/dataset/")
load("ds_hughes.rda")
load("Hughes_2009.rda")
load("hughes_2012.rda")
load("zhang_2014.rda")

## MetaCycle 
library(MetaCycle)
write.csv(hughes_2009, file="hughes_2009.csv", row.names=FALSE)
meta2d(infile="hughes_2009.csv",filestyle="csv", outdir="example",
       minper=20, maxper=28,
       timepoints=rep(seq(18, 65, by=1),each=1),
       outIntegration="onlyIntegration", ARSdefaultPer=24,
       outRawData=T)

write.csv(ds[,-2], file="ds.csv", row.names=FALSE)
meta2d(infile="ds.csv",filestyle="csv", outdir="example",
       minper=20, maxper=28,
       timepoints=rep(seq(18, 65, by=2),each=1),
       outIntegration="onlyIntegration", ARSdefaultPer=24,
       outRawData=T)

write.csv(hughes_2012[,-2], file="hughes_2012.csv", row.names=FALSE)
meta2d(infile="hughes_2012.csv",filestyle="csv", outdir="example",
       minper=20, maxper=28,
       timepoints=rep(seq(0, 46, by=2),each=1),
       outIntegration="onlyIntegration", ARSdefaultPer=24,
       outRawData=T)

write.csv(zhang[, -2], file="zhang.csv", row.names=F)
meta2d(infile="zhang.csv",filestyle="csv", outdir="example",
       minper=20, maxper=28,
       timepoints=rep(seq(18, 64, by=2),each=1),
       outIntegration="onlyIntegration", ARSdefaultPer=24,
       outRawData=T)

##BIOCYCLE 
library(readr)
write_tsv(zhang, "zhang.tsv")
write_tsv(hughes_2009, "hughes_2009.tsv")
write_tsv(hughes_2012, "hughes_2012.tsv")
write_tsv(ds, "ds.tsv")
#Rscript BioCycle.R -i /Users/mwen/Desktop/hughes_2009.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
#Rscript BioCycle.R -i /Users/mwen/Desktop/hughes_2012.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
#Rscript BioCycle.R -i /Users/mwen/Desktop/ds.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
#Rscript BioCycle.R -i /Users/mwen/Desktop/zhang.tsv -o /Users/mwen/Desktop/ -s 20  -e 28


##RAIN 
library(rain)
library("magrittr")
ds1 <- subset(ds, select=-1)
hughes12 <- subset(hughes_2012, select=c(-1,-2))
hughes09 <- subset(hughes_2009, select=c(-1,-2))
zhang1 <- subset(zhang, select=c(-1,-2))


ds_rain <- rain(t(ds[,-1]),deltat=2, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                nr.series = 1,  adjp.method = "ABH", verbose = FALSE)
hughes_2012_rain <- rain(t(hughes_2012[,-c(1,2)]),deltat=2, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                      nr.series = 1,  adjp.method = "ABH", verbose = FALSE)
hughes_2009_rain <- rain(t(hughes_2009[,-c(1,2)]),deltat=1, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                      nr.series = 1,  adjp.method = "ABH", verbose = FALSE)
zhang_rain <- rain(t(zhang[,c(1,2)]),deltat=2, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                   nr.series = 1,  adjp.method = "ABH", verbose = FALSE)

#for eJTK
write.table(ds, "ds.txt", row.names = F, col.names = T)
write.table(hughes_2009,"hughes_2009.txt", row.names = F, col.names = T)
write.table(zhang,"zhang.txt", row.names = F, col.names = T)
write.table(hughes_2012,"Hughes_2012.txt", row.names = F, col.names = T)
# ./eJTK-CalcP.py -f example/ds.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT
# ./eJTK-CalcP.py -f example/hughes_2009.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT
# ./eJTK-CalcP.py -f example/hughes_2012.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT
# ./eJTK-CalcP.py -f example/zhang.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT










