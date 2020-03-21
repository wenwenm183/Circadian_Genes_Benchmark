# read results in
meta.s8.asy.p24.t48.f8.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.asy.p24.t48.f8.new.csv')
meta.s8.nonsta.p24.t48.f8.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.nonsta.p24.t48.f8.new.csv')
meta.s8.miss1.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.miss1.p24.t48.csv')
meta.s8.miss5.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.miss8.p24.t48.csv')
meta.s8.miss10.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.miss10.p24.t48.csv')
meta.s8.snr0.5.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.snr0.5.p24.t48.new.csv')
meta.s8.snr1.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.snr1.p24.t48.new.csv')
meta.s8.snr2.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.snr2.p24.t48.new.csv')
meta.s8.snr3.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.snr3.p24.t48.new.csv')
meta.s8.uneven.p24.t48.1.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.uneven.p24.t48.1.new.csv')
meta.s8.uneven.p24.t48.2.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.uneven.p24.t48.2.new.csv')
meta.s8.uneven.p24.t48.4.new <- read.csv('./sim_12_points_new/meta_rain/meta2d_s8.uneven.p24.t48.4.new.csv')



rain.s8.asy.p24.t48.f8.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.asy.p24.t48.f8.new.csv')
rain.s8.nonsta.p24.t48.f8.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.nonsta.p24.t48.f8.new.csv')
rain.s8.miss1.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/rain.s8.miss1.p24.t48.csv')
rain.s8.miss5.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/rain.s8.miss8.p24.t48.csv')
rain.s8.miss10.p24.t48 <- read.csv('./sim_12_points_new/meta_rain/rain.s8.miss10.p24.t48.csv')
rain.s8.snr0.5.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.snr0.5.p24.t48.new.csv')
rain.s8.snr1.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.snr1.p24.t48.new.csv')
rain.s8.snr2.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.snr2.p24.t48.new.csv')
rain.s8.snr3.p24.t48.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.snr3.p24.t48.new.csv')
rain.s8.uneven.p24.t48.1.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.uneven.p24.t48.1.new.csv')
rain.s8.uneven.p24.t48.2.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.uneven.p24.t48.2.new.csv')
rain.s8.uneven.p24.t48.4.new <- read.csv('./sim_12_points_new/meta_rain/rain.s8.uneven.p24.t48.4.new.csv')


write.table(s8.nonsta.p24.t48.f8.new[,-1], './sim_12_points_new/bc/data/s8.nonsta.p24.t48.f8.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.asy.p24.t48.f8.new[,-1], './sim_12_points_new/bc/data/s8.asy.p24.t48.f8.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.uneven.p24.t48.1.new[,-1], './sim_12_points_new/bc/data/s8.uneven.p24.t48.1.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.uneven.p24.t48.2.new[,-1], './sim_12_points_new/bc/data/s8.uneven.p24.t48.2.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.uneven.p24.t48.4.new[,-1], './sim_12_points_new/bc/data/s8.uneven.p24.t48.4.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.snr0.5.p24.t48.new[,-1], './sim_12_points_new/bc/data/s8.snr0.5.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.snr1.p24.t48.new[,-1], './sim_12_points_new/bc/data/s8.snr1.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.snr2.p24.t48.new[,-1], './sim_12_points_new/bc/data/s8.snr2.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s8.snr3.p24.t48.new[,-1], './sim_12_points_new/bc/data/s8.snr3.p24.t48.new.tsv',
            sep = '\t', col.names = NA)


bc.s8.asy.p24.t48.f8.new <- read.table('./sim_12_points_new/bc/results/bc_s8.asy.p24.t48.f8.new.tsv', header = T)
bc.s8.nonsta.p24.t48.f8.new <- read.table('./sim_12_points_new/bc/results/bc_s8.nonsta.p24.t48.f8.new.tsv', header = T)
bc.s8.snr0.5.p24.t48.new <- read.table('./sim_12_points_new/bc/results/bc_s8.snr0.5.p24.t48.new.tsv', header = T)
bc.s8.snr1.p24.t48.new <- read.table('./sim_12_points_new/bc/results/bc_s8.snr1.p24.t48.new.tsv', header = T)
bc.s8.snr2.p24.t48.new <- read.table('./sim_12_points_new/bc/results/bc_s8.snr2.p24.t48.new.tsv', header = T)
bc.s8.snr3.p24.t48.new <- read.table('./sim_12_points_new/bc/results/bc_s8.snr3.p24.t48.new.tsv', header = T)
bc.s8.uneven.p24.t48.1.new <- read.table('./sim_12_points_new/bc/results/bc_s8.uneven.p24.t48.1.new.tsv', header = T)
bc.s8.uneven.p24.t48.2.new <- read.table('./sim_12_points_new/bc/results/bc_s8.uneven.p24.t48.2.new.tsv', header = T)
bc.s8.uneven.p24.t48.4.new <- read.table('./sim_12_points_new/bc/results/bc_s8.uneven.p24.t48.4.new.tsv', header = T)



write.table(s8.nonsta.p24.t48.f8.new[,-1], './sim_12_points_new/eJTK/s8.nonsta.p24.t48.f8.new.txt',
            sep = '\t', quote = F)
write.table(s8.asy.p24.t48.f8.new[,-1], './sim_12_points_new/eJTK/s8.asy.p24.t48.f8.new.txt',
            sep = '\t', quote = F)
write.table(s8.snr0.5.p24.t48.new[,-1], './sim_12_points_new/eJTK/s8.snr0.5.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s8.snr1.p24.t48.new[,-1], './sim_12_points_new/eJTK/s8.snr1.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s8.snr2.p24.t48.new[,-1], './sim_12_points_new/eJTK/s8.snr2.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s8.snr3.p24.t48.new[,-1], './sim_12_points_new/eJTK/s8.snr3.p24.t48.new.txt',
            sep = '\t', quote = F)

eJTK.s8.asy.p24.t48.f8.new <- read.table('./sim_12_points_new/eJTK/s8.asy.p24.t48.f8.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s8.nonsta.p24.t48.f8.new <- read.table('./sim_12_points_new/eJTK/s8.nonsta.p24.t48.f8.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s8.snr0.5.p24.t48.new <- read.table('./sim_12_points_new/eJTK/s8.snr0.5.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s8.snr1.p24.t48.new <- read.table('./sim_12_points_new/eJTK/s8.snr1.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s8.snr2.p24.t48.new <- read.table('./sim_12_points_new/eJTK/s8.snr2.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s8.snr3.p24.t48.new <- read.table('./sim_12_points_new/eJTK/s8.snr3.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)


auc_12_table_new <- function(){
  ds <- AUC_gp1(meta.s8.nonsta.p24.t48.f8.new, rain.s8.nonsta.p24.t48.f8.new, bc.s8.nonsta.p24.t48.f8.new, eJTK.s8.nonsta.p24.t48.f8.new, 6000, 6000)
  ds <- rbind(ds, AUC_gp1(meta.s8.asy.p24.t48.f8.new, rain.s8.asy.p24.t48.f8.new, bc.s8.asy.p24.t48.f8.new, eJTK.s8.asy.p24.t48.f8.new, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp3(meta.s8.uneven.p24.t48.1.new, rain.s8.uneven.p24.t48.1.new, bc.s8.uneven.p24.t48.1.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp3(meta.s8.uneven.p24.t48.2.new, rain.s8.uneven.p24.t48.2.new, bc.s8.uneven.p24.t48.2.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp3(meta.s8.uneven.p24.t48.4.new, rain.s8.uneven.p24.t48.4.new, bc.s8.uneven.p24.t48.4.new, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp4(meta.s8.miss1.p24.t48, rain.s8.miss1.p24.t48, 6000, 6000))
  ds <- rbind(ds, AUC_gp4(meta.s8.miss5.p24.t48, rain.s8.miss5.p24.t48, 6000, 6000))
  ds <- rbind(ds, AUC_gp4(meta.s8.miss10.p24.t48, rain.s8.miss10.p24.t48, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp1(meta.s8.snr0.5.p24.t48.new, rain.s8.snr0.5.p24.t48.new, bc.s8.snr0.5.p24.t48.new, eJTK.s8.snr0.5.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s8.snr1.p24.t48.new, rain.s8.snr1.p24.t48.new, bc.s8.snr1.p24.t48.new, eJTK.s8.snr1.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s8.snr2.p24.t48.new, rain.s8.snr2.p24.t48.new, bc.s8.snr2.p24.t48.new, eJTK.s8.snr2.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s8.snr3.p24.t48.new, rain.s8.snr3.p24.t48.new, bc.s8.snr3.p24.t48.new, eJTK.s8.snr3.p24.t48.new, 6000, 6000))
  return(ds)
}

AUC_12_table_new <- auc_12_table_new()
AUC_12_table_new$AUC <- as.numeric(as.character(AUC_12_table_new$AUC))
write.csv(AUC_12_table_new, './sim_12_points_new/AUC_12_table_new.csv')




# AUC plot
library(tidyr)
# auc.s8.group1.new <- matrix(nrow = 4, ncol = 7, dimnames = list(c('1 h/1 day', '2 h/1 day','2 h/2 days', '4 h/2 days'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
# auc.s8.group1.new[1,] <- AUC_table[c(8:14),3]
# auc.s8.group1.new[2,] <- AUC_table[c(1:7),3]
# auc.s8.group1.new[3,] <- AUC_table[c(15:21),3]
# auc.s8.group1.new[4,] <- AUC_table[c(22:28),3]
# auc.s8.group1.ds <- data.frame(auc.s8.group1.new, Category = row.names(auc.s8.group1.new),
#                                row.names = NULL)
# auc.s8.group1.ds <- gather(auc.s8.group1.ds, Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
# auc.s8.group1.ds$Method <- factor(auc.s8.group1.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))

auc.s8.group1.new <- matrix(nrow = 3, ncol = 7, dimnames = list(c('2 h/1 day','2 h/2 days', '4 h/2 days'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
auc.s8.group1.new[1,] <- AUC_table[c(1:7),3]
auc.s8.group1.new[2,] <- AUC_table[c(15:21),3]
auc.s8.group1.new[3,] <- AUC_table[c(22:28),3]
auc.s8.group1.ds <- data.frame(auc.s8.group1.new, Category = row.names(auc.s8.group1.new),
                               row.names = NULL)
auc.s8.group1.ds <- gather(auc.s8.group1.ds, Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
auc.s8.group1.ds$Method <- factor(auc.s8.group1.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))


auc.s8.group2.new <- matrix(nrow = 3, ncol = 7, dimnames = list(c('Stationary', 'Non-stationary','Asymmetric'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
auc.s8.group2.new[1,] <- AUC_table[c(1:7),3]
auc.s8.group2.new[2,] <- AUC_12_table_new[c(1:7),3]
auc.s8.group2.new[3,] <- AUC_12_table_new[c(8:14),3]
auc.s8.group2.new.ds <- gather(data.frame(auc.s8.group2.new, Category = row.names(auc.s8.group2.new),
                                      row.names = NULL), Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
auc.s8.group2.new.ds$Method <- factor(auc.s8.group2.new.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))
auc.s8.group2.new.ds$Category <- factor(auc.s8.group2.new.ds$Category, levels = c('Stationary','Non-stationary','Asymmetric'))

auc.s8.group3.new <- matrix(nrow = 3, ncol = 3, dimnames = list(c('1', '2', '3'),c('LS/MC','RAIN','BC')))
auc.s8.group3.new[1,] <- AUC_12_table_new[c(15:17),3]
auc.s8.group3.new[2,] <- AUC_12_table_new[c(18:20),3]
auc.s8.group3.new[3,] <- AUC_12_table_new[c(21:23),3]
auc.s8.group3.new.ds <- gather(data.frame(auc.s8.group3.new, Category = row.names(auc.s8.group3.new),
                                      row.names = NULL), Method, AUC,LS.MC,RAIN,BC)
auc.s8.group3.new.ds$Method <- factor(auc.s8.group3.new.ds$Method, levels = c('LS.MC','RAIN','BC'),
                                  labels = c('LS/MC','RAIN','BC'))

auc.s8.group4.new <- matrix(nrow = 3, ncol = 4, dimnames = list(c('1%', '5%', '10%'),c('LS','JTK','RAIN','MC')))
auc.s8.group4.new[1,] <- AUC_12_table_new[c(24:27),3]
auc.s8.group4.new[2,] <- AUC_12_table_new[c(28:31),3]
auc.s8.group4.new[3,] <- AUC_12_table_new[c(32:35),3]
auc.s8.group4.new.ds <- gather(data.frame(auc.s8.group4.new, Category = row.names(auc.s8.group4.new),
                                      row.names = NULL), Method, AUC,LS,JTK,RAIN,MC)
auc.s8.group4.new.ds$Method <- factor(auc.s8.group4.new.ds$Method, levels = c('LS','JTK','RAIN','MC'))
auc.s8.group4.new.ds$Category <- factor(auc.s8.group4.new.ds$Category, levels = c('1%', '5%', '10%'))

auc.s8.group5.new <- matrix(nrow = 4, ncol = 6, dimnames = list(c('2 h/1 day X 1', '4 h/1 day X 2', '2 h/2 days X 1' , '4 h/2 days X 2'),c('LS', 'JTK','RAIN','eJTK','MC','BC')))
auc.s8.group5.new[1,] <- AUC_table[c(1,3:7),3]
auc.s8.group5.new[2,] <- AUC_table[c(64:69),3]
auc.s8.group5.new[3,] <- AUC_table[c(16:21),3]
auc.s8.group5.new[4,] <- AUC_table[c(70:75),3]
auc.s8.group5.ds <- gather(data.frame(auc.s8.group5.new, Category = row.names(auc.s8.group5.new),
                                      row.names = NULL), Method, AUC, LS, JTK, RAIN,eJTK,MC,BC)
auc.s8.group5.ds$Method <- factor(auc.s8.group5.ds$Method, levels = c('LS','JTK','RAIN','eJTK','MC','BC'))


auc.s8.group6.new <- matrix(nrow = 4, ncol = 7, dimnames = list(c('SNR 0.5:1','SNR 1:1','SNR 2:1','SNR 3:1'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
auc.s8.group6.new[1,] <- AUC_12_table_new[c(36:42),3]
auc.s8.group6.new[2,] <- AUC_12_table_new[c(36:42)+7,3]
auc.s8.group6.new[3,] <- AUC_12_table_new[c(36:42)+14,3]
auc.s8.group6.new[4,] <- AUC_12_table_new[c(36:42)+21,3]
auc.s8.group6.new.ds <- gather(data.frame(auc.s8.group6.new, Category = row.names(auc.s8.group6.new),
                                      row.names = NULL), Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
auc.s8.group6.new.ds$Method <- factor(auc.s8.group6.new.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))
auc.s8.group6.new.ds$Category <- factor(auc.s8.group6.new.ds$Category, levels = c('SNR 3:1', 'SNR 2:1', 'SNR 1:1', 'SNR 0.5:1'),
                                    labels = c('3:1', '2:1', '1:1', '0.5:1'))





# plot with legend
# group 1
library(ggplot2)
pdf(file = './figures_for_paper/AUC_simulation/12 points/12points_new_with_legend.pdf', height = 3, width = 4.5)
ggplot(data = auc.s8.group1.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.65,0.35)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') 
  

# group 2
ggplot(data = auc.s8.group2.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.30)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Waveforms', shape = 'Waveforms')

# group 3
ggplot(data = auc.s8.group3.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.33)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Uneven\ntime-points', shape = 'Uneven\ntime-points')

# group 4
ggplot(data = auc.s8.group4.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.33)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'Missing values', shape = 'Missing values')

# group 5
ggplot(data = auc.s8.group5.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.35)) +
  #theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Replicates', shape = 'Replicates')

# group 6
ggplot(data = auc.s8.group6.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  # theme(legend.position = c(0.7,0.29)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  # theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'SNR', shape = 'SNR') + guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2))
dev.off()




# plot without legend
pdf(file = './figures_for_paper/AUC_simulation/12 points/12points_new_without_legend.pdf', height = 3, width = 3.5)
ggplot(data = auc.s8.group1.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())


# group 2
ggplot(data = auc.s8.group2.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Waveforms', shape = 'Waveforms') +
  theme(panel.border = element_blank())

# group 3
ggplot(data = auc.s8.group3.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Uneven\ntime-points', shape = 'Uneven\ntime-points') +
  theme(panel.border = element_blank())

# group 4
ggplot(data = auc.s8.group4.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'Missing values', shape = 'Missing values') +
  theme(panel.border = element_blank())

# group 5
ggplot(data = auc.s8.group5.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  #theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Replicates', shape = 'Replicates') +
  theme(panel.border = element_blank())

# group 6
ggplot(data = auc.s8.group6.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  # theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'SNR', shape = 'SNR') + guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2)) + 
  theme(panel.border = element_blank())
  
dev.off()








# new sampling patterns
pdf(file = './figures_for_paper/AUC_simulation/update_without_legend.pdf', height = 3, width = 3.5)
ggplot(data = auc.s4.group1.ds[auc.s4.group1.ds$Category %in% c('4 h/1 day', '8 h/2 days'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

ggplot(data = auc.s5.group1.ds[auc.s5.group1.ds$Category %in% c('3 h/1 day', '6 h/2 days'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

ggplot(data = auc.s4.group1.ds[auc.s4.group1.ds$Category %in% c('2 h/1 day', '4 h/2 days'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

ggplot(data = auc.s4.group5.ds[auc.s4.group5.ds$Category %in% c('4 h/1 day X 1', '8 h/1 day X 2'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

ggplot(data = auc.s5.group5.ds[auc.s5.group5.ds$Category %in% c('3 h/1 day X 1', '6 h/1 day X 2'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

ggplot(data = auc.s8.group5.ds[auc.s8.group5.ds$Category %in% c('2 h/1 day X 1', '4 h/1 day X 2'),], aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())

dev.off()


