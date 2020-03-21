library(ggplot2)
library(tidyverse)
library(caret)
library(mltools)
library(DescTools)

#load all the datasets -----
setwd("/Users/mwen/Desktop/Research/4_datasets/")
load('bc.ds1.rda'); ds_bio=bc.ds1.hughes_1; rm(bc.ds1.hughes_1)
load('bc.hughes2009.rda'); hughes2009_bio=bc.hughes2009; rm(bc.hughes2009)
load('bc.hughes2012.rda');hughes2012_bio=bc.hughes2012; rm(bc.hughes2012)
load('bc.norep.zhang.rda');zhang_bio=bc.norep.zhang; rm(bc.norep.zhang)
load('meta.ds1.rda');ds_meta=meta.ds1.hughes_1; rm(meta.ds1.hughes_1)
load('meta.hughes2009.rda');hughes2009_meta=meta.hughes2009; rm(meta.hughes2009)
load('meta.hughes2012.rda');hughes2012_meta=meta.hughes2012; rm(meta.hughes2012)
load('meta.norep.zhang.rda');zhang_meta=meta.norep.zhang; rm(meta.norep.zhang)
load("rain.ds1.rda");ds_rain=rain.ds1.hughes_1; rm(rain.ds1.hughes_1)
library(qvalue)
  ds_rain$qvalue=qvalue(ds_rain$pVal)$qvalue
load("rain.hughes2009.rda");hughes2009_rain= rain.hughes2009; rm(rain.hughes2009)
  hughes2009_rain$qvalue=qvalue(hughes2009_rain$pVal)$qvalue
load("rain.hughes2012.rda");hughes2012_rain =cbind(hughes2012_meta[,1],rain.hughes2012); rm(rain.hughes2012)
  hughes2012_rain$qvalue=qvalue(hughes2012_rain$pVal)$qvalue
load("rain.norep.zhang.rda");zhang_rain=rain.norep.zhang; rm(rain.norep.zhang)
  zhang_rain$qvalue=qvalue(zhang_rain$pVal)$qvalue
load("eJTK.ds1.rda"); ds_eJTK=eJTK.ds1.hughes_1; rm(eJTK.ds1.hughes_1)
load("eJTK.hughes2009.rda");hughes2009_eJTK= eJTK.hughes2009; rm(eJTK.hughes2009)
load("eJTK.hughes2012.rda");hughes2012_eJTK=eJTK.hughes2012; rm(eJTK.hughes2012)
load("eJTK.zhang.rda");zhang_eJTK=eJTK.norep.zhang; rm(eJTK.norep.zhang)



# combine all the p-value results for each dataset----
combine_dataset=function(name,...){
  name1=paste(name,"meta",sep = "_")
  name2=paste(name,"rain",sep = "_")
  name3=paste(name,"bio",sep = "_")
  name4=paste(name, "eJTK", sep="_")
  full = cbind(cbind(cbind(get(name1)[,c(1, 2, 7, 12, 17)], 
                     get(name2)[,2]),
                    get(name3)[,2]), 
                    get(name4)[,20])
  names(full) = c("ID", "ARSER", "JTK_CYCLE", "LS", "MetaCycle", "RAIN", "BIO_CYCLE", "eJTK_CYCLE")
  return(full)
}
ds_full <- combine_dataset("ds")
hughes2009_full <- combine_dataset("hughes2009")
hughes2012_full <- combine_dataset("hughes2012")
zhang_full <- combine_dataset("zhang")

#combine all the q-value results for each dataset 
qval= function(name,...) {
  name1=paste(name,"meta",sep = "_")
  name2=paste(name,"rain",sep = "_")
  name3=paste(name,"bio",sep = "_")
  name4=paste(name, "eJTK", sep="_")
  full = cbind(cbind(cbind(get(name1)[,c(1, 3, 8, 13, 18)], 
                           get(name2)[,6]),
                     get(name3)[,3]), 
               get(name4)[,21])
  names(full) = c("ID", "ARSER", "JTK_CYCLE", "LS", "MetaCycle", "RAIN", "BIO_CYCLE", "eJTK_CYCLE")
  return(full)
}

ds_q <- qval("ds")
hughes2009_q <- qval("hughes2009")
hughes2012_q <- qval("hughes2012")
zhang_q <- qval("zhang")

rm(ds_bio,ds_meta, ds_rain, ds_eJTK, hughes2009_bio, hughes2009_meta, hughes2009_rain, hughes2009_eJTK,
   hughes2012_bio, hughes2012_meta, hughes2012_rain, hughes2012_eJTK, zhang_bio, zhang_meta, zhang_rain, zhang_eJTK)

#Dark:Dark Dataset Analysis 

#Generate table of Precision and Recall----
        #import the positive and negative benchmark genes 
        library("readxl")
        circadian <- read_excel("../journal.pone.0046961.s009.xls", sheet = 1)[,1]
        noncircadian <- read_excel("../journal.pone.0046961.s009.xls", sheet = 3)
        names(circadian) <- "ID"; names(noncircadian) <- "ID"

# generate positive and negative control datasets
          # we will only be using these genes to assess ROC, AUC, precision and recall
          for (i in c("ds", "hughes2009", "hughes2012", "zhang" )) {
            a=paste(i,"full", sep = "_")
            temp_p = merge(circadian, get(a), by="ID")
            names(temp_p) <- c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
            temp_n = merge(noncircadian, get(a), by="ID")
            names(temp_n) <- c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
            assign(paste(i,"pcontrol",sep = "_"), temp_p)
            assign(paste(i,"ncontrol",sep = "_"), temp_n)
          }
          for (i in c("ds", "hughes2009", "hughes2012", "zhang" )) {
            a=paste(i,"q", sep = "_")
            temp_p = merge(circadian, get(a), by="ID")
            names(temp_p) <- c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
            temp_n = merge(noncircadian, get(a), by="ID")
            names(temp_n) <- c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
            assign(paste(a,"pcontrol",sep = "_"), temp_p)
            assign(paste(a,"ncontrol",sep = "_"), temp_n)
          }

# generate the table of measures 
    #q-value threshold
    dataset = c("hughes2009", "ds","hughes2012", "zhang")
    method = c('LS','ARS','JTK','RAIN', 'eJTK', "MC", "BC")
    table_pr = matrix(ncol = 9, nrow=12)
    colnames(table_pr) =c("Dataset", "Measures", 'LS','ARS','JTK','RAIN','eJTK', "MC", "BC")
    for (data in dataset) {
      row = which(dataset == data)
      temp_p = get(paste(data,"q_pcontrol", sep="_"))
      temp_n = get(paste(data,"q_ncontrol", sep="_"))
      for (m in method){
        col= which(method == m)
        TP <- sum(temp_p[[m]] <=0.05)
        FN <- sum(temp_p[[m]] > 0.05)
        TN <- sum(temp_n[[m]] > 0.05)
        FP <- sum(temp_n[[m]] <=0.05)
        
        table_pr[3*row, 1] =data; table_pr[3*row-1, 1] =data; table_pr[3*row -2 , 1] = data 
        table_pr[3*row, 2] ="f_measure"; table_pr[3*row-1, 2] ="recall"; table_pr[3*row -2 , 2] = "precision"
        
        Recall=TP/(TP+FN)
        Precision=TP/(TP+FP)
        table_pr[3*row, col+2]=round(2*(Recall * Precision) / (Recall + Precision), digits=4)
        table_pr[3*row-1, col+2]=round(Recall, digits = 4)
        table_pr[3*row-2, col+2]=round(Precision, digits = 4)
        
      }
    }
    prec1 = table_pr[c(1,4,7,10),] %>% data.frame() %>%
      gather(key="method", value="precision", LS:BC)
    rec1 = table_pr[c(2,5,8,11),] %>% data.frame() %>%
      gather(key="method", value="recall", LS:BC)
    fm1 = table_pr[c(3,6,9,12),] %>% data.frame() %>%
      gather(key="method", value="f_measure", LS:BC)
    q5= merge(prec1, rec1, by=c("Dataset", "method")) %>% 
      merge(fm1, by=c("Dataset", "method")) %>%
      subset( select=-c(3,5,7)) %>% 
      mutate(p_value="0.05")
    rm(prec1, rec1, fm1)

    #p-value thresholds
    threshold = c(0.000005, 0.00005, 0.0005)
    for (t in threshold) {
    for (data in dataset) {
      row = which(dataset == data)
      temp_p = get(paste(data,"pcontrol", sep="_"))
      temp_n = get(paste(data,"ncontrol", sep="_"))
      for (m in method){
      col= which(method == m)
        TP <- sum(temp_p[[m]] <=t)
        FN <- sum(temp_p[[m]] > t)
        TN <- sum(temp_n[[m]] > t)
        FP <- sum(temp_n[[m]] <=t)
        
        table_pr[3*row, 1] =data; table_pr[3*row-1, 1] =data; table_pr[3*row -2 , 1] = data 
        table_pr[3*row, 2] ="f_measure"; table_pr[3*row-1, 2] ="recall"; table_pr[3*row -2 , 2] = "precision"
        
        Recall=TP/(TP+FN)
        Precision=TP/(TP+FP)
        table_pr[3*row, col+2]=round(2*(Recall * Precision) / (Recall + Precision), digits=4)
        table_pr[3*row-1, col+2]=round(Recall, digits = 4)
        table_pr[3*row-2, col+2]=round(Precision, digits = 4)
    
      }
    }
      prec1 = table_pr[c(1,4,7,10),] %>% data.frame() %>%
        gather(key="method", value="precision", LS:BC)
      rec1 = table_pr[c(2,5,8,11),] %>% data.frame() %>%
        gather(key="method", value="recall", LS:BC)
      fm1 = table_pr[c(3,6,9,12),] %>% data.frame() %>%
        gather(key="method", value="f_measure", LS:BC)
      b1= merge(prec1, rec1, by=c("Dataset", "method")) %>% 
        merge(fm1, by=c("Dataset", "method")) %>%
        subset( select=-c(3,5,7)) %>% 
        mutate(p_value=t)
      assign(paste("b", t, sep = "_"), b1)
      rm(prec1, rec1, fm1)
    }

# precision & recall plots FIGURE 3--
library(ggplot2)

    #organize the variables for the plotting
      p_plot=rbind(`b_5e-04`, rbind(`b_5e-05`,rbind(`b_5e-06`,q5))) 
      p_plot$method <- factor(p_plot$method, levels=c('LS','ARS','JTK','RAIN', 'eJTK', "MC", "BC"))
      p_plot$Dataset = factor(p_plot$Dataset, levels=c("hughes2009", "ds", "hughes2012", "zhang"), 
                              labels = c("Hughes 2009", "Down-sampled Hughes 2009", "Hughes 2012", "Zhang 2014") )
      new_p_plot <- gather(p_plot, key="measure", value="value", precision, recall, f_measure)  %>% data.frame()
      new_p_plot$Method <- as.factor(new_p_plot$method)
      new_p_plot$measure <- factor(new_p_plot$measure, levels = c("precision", "recall","f_measure"),
                                   labels = c("Precision", "Recall", "F measure"))
      new_p_plot$p_value <- as.factor(new_p_plot$p_value)
      new_p_plot$value <- as.numeric(new_p_plot$value)
#Precision-Recall - FIGURE 3- 

pdf("../plots/p_val_plot.pdf", width = 8, height = 2)

a <- ggplot(data = new_p_plot[ which(new_p_plot$measure == "Precision"), ], aes(x = method, y = value, fill = p_value)) + 
  geom_bar(stat = 'identity',  alpha=1.0, position = position_dodge(width = 0.5)) + 
  facet_grid(~Dataset) +
  labs(x = "Method", y="Precision") +
  coord_cartesian(ylim=c(0.6,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text( color="black", 
                                   size=8, angle=90), 
        strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(face="bold")) +
  scale_fill_manual(name="p-values threshold", values=c( "#92c5de", "#5e3c99","#fdae61" ,"#d7191c"))+
  theme(axis.title.x = element_text(face= "bold", size=10))+
  theme(axis.title.y = element_text(face= "bold", size=10))+
  theme(axis.text.x= element_text(face= "bold", size=8))+
  theme(axis.text.y= element_text(face= "bold", size=8)) +
  theme(legend.title = element_text(size=8, 
                                    face="bold"))+
  theme(legend.text = element_text(size=8, 
                                   face="bold"))

b <- ggplot(data = new_p_plot[ which(new_p_plot$measure == "Recall"), ], aes(x = method, y = value, fill = p_value)) + 
  geom_bar(stat = 'identity',  alpha=1.0, position = position_dodge(width = 0.5)) + 
  facet_grid(~Dataset) +
  labs(x = "Methods", y="Recall") +
  coord_cartesian(ylim=c(0.0,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text( color="black", 
                                    size=8, angle=90), 
        strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(face="bold")) +
  scale_fill_manual(name="p-values threshold", values=c( "#92c5de", "#5e3c99","#fdae61" ,"#d7191c"))+
theme(axis.title.x = element_text(face= "bold", size=10))+
  theme(axis.title.y = element_text(face= "bold", size=10))+
  theme(axis.text.x= element_text(face= "bold", size=8))+
  theme(axis.text.y= element_text(face= "bold", size=8)) +
  theme(legend.title = element_text(size=8, 
                                    face="bold"))+
  theme(legend.text = element_text(size=8, 
                                   face="bold"))



c<- ggplot(data = new_p_plot[ which(new_p_plot$measure == "F measure"), ], aes(x = method, y = value, fill = p_value)) + 
  geom_bar(stat = 'identity',  alpha=1.0, position = position_dodge(width = 0.5)) + 
  facet_grid(~Dataset) +
  labs(x = "Methods", y="F measure") +
  coord_cartesian(ylim=c(0.0,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text( color="black", 
                                    size=8, angle=90), 
        strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(face="bold")) +
  scale_fill_manual(name="p-values threshold", values=c( "#92c5de", "#5e3c99","#fdae61" ,"#d7191c"))+
  theme(axis.title.x = element_text(face= "bold", size=10))+
  theme(axis.title.y = element_text(face= "bold", size=10))+
  theme(axis.text.x= element_text(face= "bold", size=8))+
  theme(axis.text.y= element_text(face= "bold", size=8)) +
  theme(legend.title = element_text(size=8, 
                                    face="bold"))+
  theme(legend.text = element_text(size=8, 
                                   face="bold"))
print(a)
print(b)
print(c)
dev.off() 

rm(`b_5e-04`, `b_5e-05`, `b_5e-06`, b1, new_p_plot, p_plot, q5, table_pr, table_pr1, temp_n, temp_p)


#ROC Curves - FIGURE 2- ----- 
#generate table of sensitivity and 1- specificity 

dataset = c("hughes2009", "ds","hughes2012", "zhang")
method = c('LS','ARS','JTK','RAIN', 'eJTK', "MC", "BC")
table_roc = matrix(ncol = 14, nrow=1001)
colnames(table_roc) = c( "LS", "LS_sp", "ARS","ARS_sp",  "JTK","JTK_sp","RAIN", "RAIN_sp", 
                         "eJTK", "eJTK_sp", "MC","MC_sp", "BC", "BC_sp")
for (data in dataset) {
  temp_p = get(paste(data,"pcontrol", sep="_"))[, -1]
  temp_n = get(paste(data,"ncontrol", sep="_"))[, -1]
  
  for (i in method){
    col=which(method == i)
    x = seq(0.0, 1.00, by=0.001)
    for (k in 1:1001){
      # how we define TP <= or just = 
    TP <- nrow(subset(temp_p[temp_p[[i]] <= x[k], ], select=1))
    FN <- nrow(subset(temp_p[temp_p[[i]] > x[k], ], select=1))
    FP <- nrow(subset(temp_n[temp_n[[i]] <= x[k], ], select=1))
    TN <- nrow(subset(temp_n[temp_n[[i]] > x[k], ], select=1))
    # sensitivity 
    table_roc[k,2*col-1] =TP/(TP+FN)
    #1 - specificity
    table_roc[k,2*col] =1-TN/(TN+FP)
    
    }
  }
  
  temp = matrix(data=0, nrow=1, ncol=14)
  colnames(temp)= c( "LS", "LS_sp", "ARS","ARS_sp",  "JTK","JTK_sp","RAIN", "RAIN_sp", 
                     "eJTK", "eJTK_sp", "MC","MC_sp", "BC", "BC_sp")
  table_roc=rbind(temp, table_roc)
  
  sensitivity = table_roc %>% data.frame() %>%
    subset(select=c(seq(1,14,by=2))) %>% 
    gather(key="Method", value="sensitivity", LS:BC)
  
  specificity= table_roc %>% data.frame() %>%
    subset(select=c(seq(2,14,by=2))) %>% 
    gather(key="Method_1", value="specificity", LS_sp:BC_sp)
  
  temp = cbind(sensitivity, specificity)
  assign(paste(data, "roc",sep="_"), temp)
  
}

#plot ROC curve 
a1= matrix(ncol=7, nrow = 1)
for (data in dataset) {
  temp = get(paste(data, "roc", sep = "_"))
  setwd("../plots/")
  filename= paste(data, "pdf", sep = ".")
  method = c('LS','ARS','JTK','RAIN','eJTK', "MC", "BC")
  #label with AUC
  for (m in method) {
    j= which(method == m)
    auc=filter(temp, Method == m)
    library(DescTools)
      a=round(AUC(auc[,4], auc[,2]), 4)
      a1[j] =paste(m,"(AUC:",a,")")
  }
  
  temp$Method1 = factor(temp$Method, levels=c('LS','ARS','JTK','RAIN','eJTK', "MC", "BC"))
  #make a pdf file
  pdf(filename, width=3.25, height=3.25)
 p <- ggplot(temp, mapping=aes(x=specificity, y=sensitivity)) + 
    geom_line(aes(colour=Method1)) + 
    labs(x='1 - Specificity',
         y='Sensitivity')+ 
   scale_colour_discrete(name = "Methods", labels = a1) +
   theme_bw() +
   theme(plot.background = element_blank()
         ,panel.grid.major = element_blank()
         ,panel.grid.minor = element_blank()) +
   theme(plot.title = element_text(hjust = 0.5)) +
   theme(plot.title = element_text(size=20))+ 
   coord_fixed()+ 
   theme(axis.title.x = element_text(face= "bold", size=10))+
   theme(axis.title.y = element_text(face= "bold", size=10))+
   theme(axis.text.x= element_text(face= "bold", size=8))+
   theme(axis.text.y= element_text(face= "bold", size=8)) +
   theme(legend.title = element_text(size=7, 
                                     face="bold"))+
   theme(legend.text = element_text(size=6, 
                                    face="bold")) +
   theme(legend.position = c(0.75, 0.28)) + 
   theme(legend.key.height=unit(0.8,"line"))+ 
   theme(legend.key.width=unit(1,"line"))
   
    print(p)
  dev.off() 
}

rm(ds_roc, hughes2012_roc, hughes2009_roc, zhang_roc, table_roc, sensitivity, specificity, temp, temp_n, temp_p)









rm(ds_ncontrol, ds_pcontrol, hughes2009_ncontrol, hughes2009_pcontrol, hughes2012_pcontrol, 
   hughes2012_ncontrol, zhang_ncontrol, zhang_pcontrol)
rm(a1, auc, b, c, ds_q, ds_q_ncontrol, ds_q_pcontrol, hughes2009_q, hughes2009_q_ncontrol, 
   hughes2009_q_pcontrol,hughes2012_q, hughes2012_q_ncontrol, hughes2012_q_pcontrol, zhang_q, 
   zhang_q_pcontrol, zhang_q_ncontrol,p)
#Venn FIGURE 3A -----
intersect <- intersect(intersect(intersect(ds_full$ID, hughes2009_full$ID), hughes2012_full$ID), zhang_full$ID) %>% data.frame()
names(intersect) ="ID"

  names(ds_full)= c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  names(hughes2009_full)= c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  names(hughes2012_full)= c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  names(zhang_full)= c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
  dataset = c("hughes2009", "ds","hughes2012", "zhang")
  for (d in dataset) {
    library(qvalue)
    data=get(paste(d,"full", sep = "_"))
    temp=merge(intersect, data, by="ID")
    temp1=temp[,1] %>% data.frame()
    for (num in c(2:8)) {
      temp1[,num]= qvalue(p=temp[,num])$qvalue
    }
    names(temp1) <- c("ID", "ARS", "JTK", "LS", "MC", "RAIN", "BC", "eJTK")
    assign(paste("q", d, sep = "_"), temp1)
  }


  method = c('LS','ARS','JTK','RAIN', 'eJTK', "MC", "BC")
  for (m in method) {
    setwd("/Users/mwen/Desktop/plots/")
    library(VennDiagram)
    x1 = subset(q_hughes2009[q_hughes2009[[m]] <= 0.05,], select=1)
    x2 = subset(q_hughes2012[q_hughes2012[[m]] <= 0.05,], select=1)
    x3 = subset(q_ds[q_ds[[m]] <= 0.05,], select=1)
    x4 = subset(q_zhang[q_zhang[[m]] <= 0.05,], select=1)
    
    venn.plot <- venn.diagram(
      x = list(
        Hughes = x1$ID,
        Hughes_2012= x2$ID,
        Downsampled= x3$ID,
        Zhang= x4$ID
      ),
      filename = NULL,
      main = m,
      main.cex = 4,
      col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
      category = c("Hughes 2009", "Hughes 2012", "Downsampled Hughes 2009", "Zhang 2014"),
      cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
      cex = c(1, 1, 1, 1,  1, 2.5, 1, 1, 1, 1,  1, 1, 1, 1, 1),
      lty =rep("solid"),
      cat.fontfamily = "serif",
      cat.cex = 1.0,
      cat.fontface = "bold",
      margin = 0.1,
      ext.text= TRUE
    )
    

    library(grDevices)
    title=paste(m, "pdf", sep = ".")
    pdf(file=title)
    grid.draw(venn.plot)
    dev.off()
  
  }


#JACCARD and SORENSEN 3B -----
  
  library(OmicsMarkeR)
  for (m in method) {
    setwd("/Users/mwen/Desktop/plots/")

    x1 = subset(q_hughes2009[q_hughes2009[[m]] <= 0.05,], select=1) %>%
      unlist(use.names=FALSE) %>% as.character()
    x3 = subset(q_hughes2012[q_hughes2012[[m]] <= 0.05,], select=1)%>%
      unlist(use.names=FALSE) %>% as.character()
    x2 = subset(q_ds[q_ds[[m]] <= 0.05,], select=1) %>%
      unlist(use.names=FALSE) %>% as.character()
    x4 = subset(q_zhang[q_zhang[[m]] <= 0.05,], select=1) %>%
      unlist(use.names=FALSE) %>% as.character()
    
    #plot
    library(RColorBrewer)
    library(pheatmap)
    library(grid)
    J=matrix(nrow=4, ncol=4, data=0)
    colnames(J)=rownames(J)=c('Hughes2009','Downsample','Hughes2012','Zhang2014')
    diag(J)=1
    #calculate the jaccard index value
    J[1,2:4]=c(round(jaccard(x1, x2),digits=4), round(jaccard(x1, x3),digits=4),round(jaccard(x1, x4),digits=4))
    J[2,3:4]=c(round(jaccard(x2, x3),digits=4),round(jaccard(x2, x4),digits=4))
    J[3,4]=  round(jaccard(x3, x4),digits=4)
    name=paste(paste(m, "jaccard", sep="-"), "pdf", sep=".")
    pdf(file=name, width=3, height=3)
    pheatmap(J, cluster_rows = F, cluster_cols = F,
             color = colorRampPalette(rev(c("red","white")))(50),
             show_rownames = F, show_colnames = F, legend = F)
    grid.text(J[1,2:4], x=0.125, y=rev(seq(0,0.5,0.25)+0.125), gp=gpar(fontsize=8))
    grid.text(J[2,3:4], x=0.25+0.125, y=rev(seq(0,0.25,0.25)+0.125), gp=gpar(fontsize=8))
    grid.text(J[3,4], x=0.5+0.125, y=0.125, gp=gpar(fontsize=8))
    grid.text(c('Hughes2009','Downsample','Hughes2012','Zhang2014'), x=seq(0,0.75,0.25)+0.125, y=rev(seq(0,0.75,0.25)+0.125), gp=gpar(fontsize=8))
    dev.off()
    
  }
  
  
  # Sorensen index 
  
  method = c('LS','ARS','JTK','RAIN', 'eJTK', "MC", "BC")
  soren1=matrix(data=NA, nrow=7, ncol=2)
  for (m in method) {
    row = which(method == m)
    x1 = subset(q_hughes2009[q_hughes2009[[m]] <= 0.05,], select=1)
    x2 = subset(q_ds[q_ds[[m]] <= 0.05,], select=1)
    x3 = subset(q_hughes2012[q_hughes2012[[m]] <= 0.05,], select=1)
    x4 = subset(q_zhang[q_zhang[[m]] <= 0.05,], select=1)
    
    AB <- intersect(x1,x2)
    AC <- intersect(x1,x3)
    AD <- intersect(x1,x4)
    BC <- intersect(x2,x3)
    BD <- intersect(x2,x4)
    CD <- intersect(x3,x4)
    ABC<- intersect(intersect(x1,x2),x3)
    BCD <- intersect(intersect(x2,x3), x4)
    ACD <- intersect(intersect(x1,x3),x4)
    ABD <- intersect(intersect(x1,x2),x4)
    ABCD <- intersect(intersect(intersect(x1,x2),x3), x4)
    
    y <- nrow(AB)+nrow(BC)+nrow(AD)+nrow(AC)+nrow(BD)+nrow(CD)-nrow(ABC)-nrow(BCD)-nrow(ACD)-nrow(ABD)+nrow(ABCD)
    z <- nrow(x1)+nrow(x2)+nrow(x3)+nrow(x4)
    
    soren1[row, 1]=m
    soren1[row, 2]=round(y/z*(4/3), digits=4)
    
  }
  
  rm(x1, x2, x3, x4,AB, ABC, ABCD, AC, AD, BC, BD, BCD, CD,ABD, ACD,intersect,venn.plot,
     J, temp1, q_hughes2009, q_hughes2012, q_ds, q_zhang,temp)
  
#FLAW WITH FDR -----
setwd("/Users/mwen/Desktop/Research/dataset/")
load("Hughes_2009.rda")
hughes <- hughes_2009[,-c(1,2)]
#permutate
  library(vegan)
  library(tidyverse)
    temp <- permatfull(hughes, times=1, fixedmar = "rows", shuffle = "samp")
    temp1 <- temp[["perm"]][[1]] %>% data.frame()
    hughes_shuffle <- cbind(hughes_2009$Gene.Symbol, temp1)
    rm(temp, temp1)
    
# RUN WITH DIFFERENT METHODS ---
  #MetaCycle
    setwd("/Users/mwen/Desktop/Research/shuffle/")
    library(readr)
    write.csv(hughes_shuffle, file="hughes_shuffle.csv", row.names=F)
    library(MetaCycle)
    meta2d(infile="hughes_shuffle.csv",filestyle="csv",outdir="MetaCycle",
           minper=20, maxper=28,
           timepoints=seq(18, 65, by=1),
           outIntegration="onlyIntegration", ARSdefaultPer=24,
           outRawData=FALSE)
  #RAIN 
    library(rain)
    library("magrittr")
    hughes_2009_rain <- rain(t(hughes_shuffle[,-c(1)]),deltat=1, period=24, period.delta = 4, peak.border = c(0.3, 0.7),
                             nr.series = 1,  adjp.method = "ABH", verbose = FALSE)
    #save(hughes_2009_rain, file="hughes_2009_rain.rda")

  #BC 
    write_tsv(hughes_shuffle, "hughes_shuffle.tsv")
    #Rscript BioCycle.R -i /Users/mwen/Desktop/hughes_shuffle.tsv -o /Users/mwen/Desktop/ -s 20  -e 28
    
  #eJTK 
    write.table(hughes_shuffle, "hughes_shuffle.txt", row.names = F, col.names = T)
    # ./eJTK-CalcP.py -f example/hughes_shuffle.txt -w ref_files/waveform_cosine.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2_OTHERTEXT
    
    
  #Results with shuffled --
  setwd("/Users/mwen/Desktop/Research/shuffle/")
    mc <- read.csv("meta2d_hughes_shuffle.csv", header=TRUE, sep=",") 
    eJTK <- read.table(file="hughes_shuffle_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T, sep = "\t") 
    bc =read_tsv(file="hughes_biocycle.tsv", col_names=T) 
    load("hughes_2009_rain.rda"); rain=hughes_2009_rain; rm(hughes_2009_rain)
    # histograms for the shuffled results 
  pval=cbind(cbind(cbind(mc[,c(2,7,12,17)],bc[,2]),rain[,1]),eJTK[,19])
  names(pval) = c("ARS", "JTK", "LS", "MC", "BC", "RAIN", "eJTK")
    method = c("ARS", "JTK", "LS", "MC", "BC", "RAIN", "eJTK")
  for(m in method) {
    setwd("/Users/mwen/Desktop/plots/")
    name=paste(paste(m,"shuffled", sep="-"),"pdf",sep=".")
    methods=which(method == m)
    pdf(name)
    hist(pval[,methods], xlab="p-values", main=m)
    dev.off()
    
  }
    
    #histograms for the original p-values of Hughes 2009 
    hughes_pval=cbind(cbind(cbind(hughes2009_meta[,c(2,7,12,17)],hughes2009_bio[,2]),hughes2009_rain[,2]),hughes2009_eJTK[,19])
    names(hughes_pval) = c("ARS", "JTK", "LS", "MC", "BC", "RAIN", "eJTK")
    method = c("ARS", "JTK", "LS", "MC", "BC", "RAIN", "eJTK")
    
    for(m in method) {
      setwd("/Users/mwen/Desktop/plots/")
      name=paste(paste(m,"original", sep="-"),"pdf",sep=".")
      methods=which(method == m)
      pdf(name)
      hist(hughes_pval[,methods], xlab="p-values", main=m)
      dev.off()
    }
    
    # before and after plot of Cry1
    hughes_shuffle = read.csv("hughes_shuffle.csv", header=TRUE, sep=",") 
    name=c("id", "key", "value","group")
    temp= hughes_shuffle[which(hughes_shuffle$ID =="Cry1"),] %>%
      gather("key", "value", CT.18:CT.65) 
    temp$group="after"
    names(temp)=name
    temp1 <- hughes_2009[which(hughes_2009$Gene.Symbol=="Cry1"), -1] %>%
      gather("key", "value", GSM301348:GSM301395) 
    temp1$group="before"
    names(temp1)=name
    
    plot = rbind(temp, temp1) 
    plot$time=rep(seq(18,65,1), time=2)
    plot$group=factor(plot$group, levels=c("before", "after"))
    
    setwd("/Users/mwen/Desktop/plots/")
    pdf("before_after_shuffle.pdf", width = 4.5, height = 2.5)
    
    p <- ggplot(plot[which(plot$group=="before"),], aes(time, value))+
      geom_point(color="#e41a1c") +
      labs(x='CT (Hours)',
           y='Expression Level', 
           title = "Gene Cry1 Before Time Shuffle") +
      theme_bw() +
      theme(plot.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()) +
      theme(legend.position = "none")+
      theme(plot.title=element_text(size=12, face="bold.italic")) + 
      theme(axis.title.x = element_text(face= "bold", size=10))+
      theme(axis.title.y = element_text(face= "bold", size=10))+
      theme(axis.text.x= element_text(face= "bold", size=8))+
      theme(axis.text.y= element_text(face= "bold", size=8)) 
    
    q <- ggplot(plot[which(plot$group=="after"),], aes(time, value))+
      geom_point(color='#377eb8') +
      labs(x='CT (Hours)',
           y='Expression Level', 
           title = "Gene Cry1 After Time Shuffle") +
      theme_bw() +
      theme(plot.background = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()) +
      theme(legend.position = "none")+
      theme(plot.title=element_text(size=12, face="bold.italic")) + 
      theme(axis.title.x = element_text(face= "bold", size=10))+
      theme(axis.title.y = element_text(face= "bold", size=10))+
      theme(axis.text.x= element_text(face= "bold", size=8))+
      theme(axis.text.y= element_text(face= "bold", size=8)) 
    
    print(p)
    print(q)
    dev.off() 
    

    
    
    #HEATMAP Plot 
    #BEFORE SHUFFLING 
    Hughes_2009 <- read.csv('../data/Hughes_2009.csv')
    sort(sapply(t(Hughes_2009[,-1]),sd), decreasing = T)
    a <- apply(t(Hughes_2009[,-1]),2,sd)
    order(a, decreasing = T)[1:2000]
    corr1 <- cor(t(Hughes_2009[order(a, decreasing = T)[1:2000],-1]))
    pdf(file = 'htmp.pdf', width = 2, height = 2)
    pheatmap(corr1[c(1:200),c(1:200)], show_colnames = FALSE,
             show_rownames = FALSE, cluster_rows = F, cluster_cols = F)
    dev.off()
    
    #AFTER SHUFFLING 
    hughes_shuffle <- read.csv('../data/hughes_shuffle.csv')
    sort(sapply(t(hughes_shuffle[,-1]),sd), decreasing = T)
    a <- apply(t(hughes_shuffle[,-1]),2,sd)
    order(a, decreasing = T)[1:2000]
    cor_shuffle <- cor(t(hughes_shuffle[order(a, decreasing = T)[1:2000],-1]))
    pdf(file = 'htmp_shuffled.pdf', width = 2, height = 2)
    pheatmap(cor_shuffle[c(1:200),c(1:200)], show_colnames = FALSE,
             show_rownames = FALSE, cluster_rows = F, cluster_cols = F)
    dev.off()
  

    
# PLOT GENES FOR FIGURE 1-----
    load('hughes2009.rda')
    hughes2009=hughes.new; rm(hughes.new)
    load('hughes2012.rda')
    hughes2012=hughes.2012.wenwen; rm(hughes.2012.wenwen)
    hughes2012=hughes2012[,-2]
    colnames(hughes2012)[1]='symbol'
    load('zhang2014.rda')
    zhang=norep.zhang; rm(norep.zhang)
    
    # take gene intersections
    genes=intersect(intersect(hughes.new$symbol, hughes.2012$Gene.Symbol), norep.zhang$symbol)
    
    hughes2009=as.matrix(hughes.new[match(genes, hughes.new$symbol),2:ncol(hughes.new)])
    hughes2012=as.matrix(hughes.2012[match(genes, hughes.2012$Gene.Symbol),3:ncol(hughes.2012)])
    zhang1=as.matrix(norep.zhang[match(genes, norep.zhang$symbol),2:ncol(norep.zhang)])
    
    gene='Polrf'
    gene='Clock'
    
    
    library(ggplot2)
    plot.gene=function(gene,...){
      temp1=scale(hughes2009[which(genes==gene),])
      temp2=scale(hughes2012[which(genes==gene),])
      temp3=scale(zhang1[which(genes==gene),])
      
      df=data.frame(CT=c(seq(18,65,1), seq(0,46,2),seq(18,64,2)),
                    exp=c(temp1, temp2, temp3),
                    Study=c(rep('Hughes2009',length(temp1)), rep('Hughes2012',length(temp2)), rep('Zhang2014',length(temp3))))
      ggplot(df, aes(x=CT, y=exp, group=Study)) +
        geom_line(aes(color=Study))+
        geom_point(aes(shape=Study, col=Study))+ 
        ylab('Scaled expression')+
        labs(title=c(gene))+
        theme_bw() + theme(panel.border = element_blank()) +
        theme(plot.title=element_text(face='bold.italic',size = 12)) + 
        theme(axis.title.x=element_text(face='bold',size = 10)) + 
        theme(axis.text.x=element_text(face='bold',size = 10)) + 
        theme(axis.title.y=element_text(face='bold',size = 10)) + 
        theme(axis.text.y=element_text(face='bold',size = 10)) + 
        theme(legend.text=element_text(face='bold',size = 10)) + 
        theme(legend.title=element_text(face='bold',size = 10)) + 
        ggsave(paste(gene,'.pdf',sep=''), width = 4, height=2)
    }
    
    
    genes[which(genes %in% gene.symbl.n$ID)]
    
    
    
    plot.gene('Mtf1')
    plot.gene('Nub1')
    plot.gene('Mynn')
    plot.gene('Snx12')
    
    plot.gene('Clock')
    plot.gene('Cry1')
    plot.gene('Npas2')
    plot.gene('Per1')
    
    plot.gene('Utp6')
    plot.gene('Mtf1')
    plot.gene('Cln3')
    plot.gene('Abcd4')
    
    
    plot.gene('Dbp')
    plot.gene('Hspa5')
    
  
  
  
  
  
  
  

  
  
  
  
  