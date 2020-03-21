library(ggplot2)
library(tidyverse)
library(caret)
library(mltools)
library(DescTools)

setwd("/Users/mwen/Desktop/Research/LD/")
load('XR.rda')
load('XR_rain.rda')
load('XR_metacycle.rda')
load('XR_eJTK.rda')
load('XR_biocycle.rda')
load('RNA.rda')
load('RNA_rain.rda')
load('RNA_metacycle.rda')
load('RNA_eJTK.rda')
load('RNA_biocycle.rda')
load('Nascent.rda')
load('Nascent_rain.rda');nascent_rain=Nascent_rain; rm(Nascent_rain)
load('Nascent_metacycle.rda')
load('Nascent_eJTK.rda')
load('Nascent_biocycle.rda')
load('gro.rda')
load('gro_rain.rda')
load('gro_metacycle.rda')
load('gro_eJTK.rda')
load('gro_biocycle.rda')


combine_dataset=function(name,...){
  name1=paste(name,"meta2d",sep = "_")
  name2=paste(name,"rain",sep = "_")
  name3=paste(name,"bio",sep = "_")
  name4=paste(name, "eJTK", sep="_")
  full = cbind(cbind(cbind(get(name1)[,c(1, 2, 7, 12)], 
                           get(name2)[,1]),
                     get(name3)[,2]), 
               get(name4)[,19])
  names(full) = c("ID", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  return(full)
}

nascent_full <- combine_dataset("nascent")
RNA_full <- combine_dataset("RNA")
XR_full <- combine_dataset("XR")
gro_full <- cbind(cbind(cbind(gro_meta2d[,c(1, 7, 12,17)], gro_rain[,1]),gro_bio[,2]),gro_eJTK[,19])
names(gro_full) = c("ID", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")

rm(gro_bio, gro_eJTK,gro_meta2d, gro_rain, gro1, nascent_bio, nascent_eJTK, nascent_meta2d, 
   nascent_rain, Nascent1, RNA_bio, RNA_eJTK, RNA_meta2d, RNA_rain, RNA1, XR, XR_bio, XR_eJTK, 
   XR_rain, XR_meta2d)

library("readxl")
circadian <- read_excel("../journal.pone.0046961.s009.xls", sheet = 1)[,1]
names(circadian) <- "ID"

for (i in c("gro", "nascent", "RNA", "XR" )) {
  a=paste(i,"full", sep = "_")
  temp_p = merge(circadian, get(a), by="ID")
  names(temp_p) <- c("ID", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  assign(paste(i,"pcontrol",sep = "_"), temp_p)
}

Nascent_pcontrol = nascent_pcontrol
GRO_pcontrol=gro_pcontrol
plot.bee=function(data,...){
  library(ggbeeswarm)
  a=paste(data,"pcontrol", sep="_")
  x=get(a)
  temp = gather(x, key="Method", value="p_Val", JTK:eJTK) 
  temp$logP = -log10(temp$p_Val)
  temp$Method <- factor(temp$Method, levels=c('LS','JTK','RAIN', 'eJTK', "MC", "BC"))
 
  p <- ggplot(temp,aes(Method, logP, color=Method)) + 
    scale_y_continuous(limits = c(0, 6)) +
    geom_quasirandom(size=.5) +
    labs(x='Method',
         y='-Log P-values',
         title= paste(paste(paste(data, "seq", sep="-"),nrow(x),sep=" "),"Positive Controls",sep=" ")) +
    theme_minimal()+
    theme(legend.position = "none", plot.title =element_text(size=14, face="bold.italic",hjust = 0.5)) +
    theme(axis.title.x = element_text(face= "bold", size=12))+
    theme(axis.title.y = element_text(face= "bold", size=12))+
    theme(axis.text.x= element_text(face= "bold", size=10))+
    theme(axis.text.y= element_text(face= "bold", size=10))
   

  return(p)
  dev.off()
  
}

setwd("/Users/mwen/Desktop/")
pdf("beeplot.pdf", width = 6, height = 4)
plot.bee("XR")
plot.bee("Nascent")
plot.bee("RNA")
plot.bee("GRO")


dev.off() 


###### of significant genes ------

setwd("/Users/mwen/Desktop/Research/LD/full/")
load("gro_full.rda"); gro=gro_full; rm(gro_full)
load("nascent_full.rda"); nascent=nascent_full; rm(nascent_full)
load("RNA_full.rda"); rna=RNA_full; rm(RNA_full)
load("XR_full.rda"); xr=XR_full; rm(XR_full)

genes=intersect(intersect(intersect(gro$ID, nascent$ID), xr$ID), rna$ID)
length(genes)

# interested genes
gro=gro[is.element(gro$ID, genes),]
nascent=nascent[is.element(nascent$ID, genes),] 
rna=rna[is.element(rna$ID, genes),]
xr=xr[is.element(xr$ID, genes),]

gro=gro[, c(1,2, 4, 8, 10, 12, 14)]
nascent=nascent[, c(1,seq(2,12,by=2))]
rna=rna[, c(1,seq(2,12,by=2))]
xr=xr[, c(1,seq(2,12,by=2))]

names(gro) <- c("ID", "BC", "RAIN", "JTK", "LS", "MC", "eJTK")
names(rna) <- c("ID", "BC", "RAIN", "JTK", "LS", "MC", "eJTK")
names(xr) <- c("ID", "BC", "RAIN", "JTK", "LS", "MC", "eJTK")
names(nascent) <- c("ID", "BC", "RAIN", "JTK", "LS", "MC", "eJTK")

dataset = c("gro", "nascent","rna", "xr")
method = c('LS','JTK','RAIN', 'eJTK', "MC", "BC")
for (data in dataset) {
  library(qvalue)
  temp=get(data)
  q = matrix(ncol=7, nrow = 9481) %>% data.frame()
  colnames(q)=c( "ID", 'LS','JTK','RAIN', 'eJTK', "MC", "BC")
  q[,1]=temp[,1]
  for (m in method) {
    col=which(method == m)
    q[,col+1]=qvalue(temp[[m]])$qvalues 
  }
  assign(paste(data, "q",sep="_"), q)
  
}


dataset = c("gro", "nascent","rna", "xr")
sig = matrix(ncol=5, nrow = 6) %>% data.frame()
colnames(sig)=c( "Method", 'GRO-seq','Nascent-seq','RNA-seq','XR-seq')
for (data in dataset) {
  temp=get(paste(data, "q",sep="_"))
  col=which(dataset == data)
  for (m in method) {
    row=which(method == m)
    sig[row,1]=m
    sig[row,col+1]=sum(temp[[m]]<=0.05)
  }
}

library(gridExtra)
pdf('sig.pdf')
grid.table(sig)
dev.off()

library("readxl")
genelist <- read.csv(file = "GeneID_symbol_convert_v1.csv",na.strings=c("","NA"), header = F,
                     stringsAsFactors = F)

sig_bio.genes=genes[gro_bio$Q_VALUE<=0.05 & Nascent_bio$Q_VALUE<=0.05 & XR_bio$Q_VALUE<=0.05 & RNA_bio$Q_VALUE<=0.05]

#unfinished because different length of the columns
dataset = c("gro", "nascent","rna", "xr")
sig_list = matrix(ncol=10, nrow = 2000) %>% data.frame()
for (data in dataset) {
  temp = get(paste(data,"q", sep = "_"))
for (m in method) {
  a =temp[temp[[m]]<=0.05,1] 
  b=genelist[!is.na(match(genelist[,2], a)),1]
  c=paste(paste(m, "ID",sep="_"),'txt', sep=".")
  
}
}




subset(xr_q, JTK < 0.05)$ID
subset(genelist, V2 = subset(xr_q, JTK < 0.05)$ID)
View(genelist[genelist$V2 %in% subset(xr_q, JTK < 0.05)$ID, ])

write.table(b, file=c, quote = F, row.names = F, col.names = F)

for (m in method) {
  a =rna_q[rna_q[[m]]<=0.05,1] 
  b=genelist[!is.na(match(genelist[,2], a)),1]
  c=paste(paste(m, "ID",sep="_"),'txt', sep=".")
  write.table(b, file=c, quote = F, row.names = F, col.names = F)
}

for (m in method) {
  a =nascent_q[nascent_q[[m]]<=0.05,1] 
  b=genelist[!is.na(match(genelist[,2], a)),1]
  c=paste(paste(m, "ID",sep="_"),'txt', sep=".")
  write.table(b, file=c, quote = F, row.names = F, col.names = F)
}

for (m in method) {
  a =gro_q[gro_q[[m]]<=0.05,1] 
  b=genelist[!is.na(match(genelist[,2], a)),1]
  c=paste(paste(m, "ID",sep="_"),'txt', sep=".")
  write.table(b, file=c, quote = F, row.names = F, col.names = F)
}


intersect_en=genelist[!is.na(match(genelist[,2], genes)),1] %>% data.frame()
write.table(intersect_en, file="inter_ENS.txt", quote = F, row.names = F, col.names = F)

which(duplicated(genelist[,2]), arr.ind = TRUE)
genelist=genelist[!is.na(genelist[,2]),c(1,2)]


#####supplements -----
library(gridExtra)
setwd("/Users/mwen/Desktop/sup/")
gro_MC <- read.delim("gro_MC.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>% 
  mutate(Method="GRO-seq MC") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>%
  select(Method, everything())
nascent_eJTK <- read.delim("nascent_eJTK.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>% 
  mutate(Method="Nascent-seq eJTK") %>% subset(Benjamini <= 0.05, select=c(1:7))   %>% select(Method, everything())
rna_BC <- read.delim("rna_BC.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>% 
  mutate(Method="RNA-seq BC") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())
rna_eJTK<- read.delim("rna_eJTK.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>%
  mutate(Method="RNA-seq eJTK") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())
rna_JTK <- read.delim("rna_JTK.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>%
  mutate(Method="RNA-seq JTK") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())
rna_MC <- read.delim("rna_MC.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>%
  mutate(Method="RNA-seq MC") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())
xr_BC <- read.delim("xr_BC.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>%
  mutate(Method="XR-seq BC") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())
xr_eJTK <- read.delim("xr_eJTK.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)] %>%
  mutate(Method="XR-seq eJTK") %>% subset(Benjamini <= 0.05, select=c(1:7)) %>% select(Method, everything())

supp_table <- rbind(rbind(rbind(rbind(rbind(rbind(rbind(gro_MC, nascent_eJTK), rna_JTK), rna_eJTK), rna_MC), rna_BC),xr_eJTK),xr_BC)
names(supp_table) <- c("Method", "Category", "Term", "Count", "Percent", "P-Value", "Benjamini" )

pdf("supp.table.pdf", height = 25, width = 15)
grid.table(supp_table)
dev.off()


dataset = c("nascent")
method=c("eJTK")
#method =c("eJTK", "BC","MC", "JTK")
for (data in dataset) {
  for (m in method) {
    file1 = paste(paste(data,m, sep = "_"),"txt", sep = ".")
    name1= paste(paste(data,m, sep = "_"),"pdf", sep = ".")
    temp <- read.delim(file1, header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)]
    temp= temp[which(temp[,6] <=0.05),]
    names(temp) <- c("Category", "Term", "Count", "Percent", "P-Value", "Benjamini" )
    pdf(name1, height = 25, width = 15)
    grid.table(temp)
    dev.off()
    }
}

pdf('rna_eJTK.pdf', height = 25, width = 15)
rna_eJTK <- read.delim("RNA_eJTK.txt", header = TRUE, sep = "\t",dec = ".")[, c(1,2,3,4,5,12)]
names(rna_eJTK) <- c("Category", "Term", "Count", "Percent", "P-Value", "Benjamini" )
rna_eJTK= rna_eJTK[which(rna_eJTK$Benjamini <=0.05),]
grid.table(rna_eJTK)
dev.off()








