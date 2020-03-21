# read results in
meta.s6.asy.p24.t48.f8.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.asy.p24.t48.f8.new.csv')
meta.s6.nonsta.p24.t48.f8.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.nonsta.p24.t48.f8.new.csv')
meta.s6.miss1.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.miss1.p24.t48.csv')
meta.s6.miss5.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.miss5.p24.t48.csv')
meta.s6.miss10.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.miss10.p24.t48.csv')
meta.s6.snr0.5.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.snr0.5.p24.t48.new.csv')
meta.s6.snr1.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.snr1.p24.t48.new.csv')
meta.s6.snr2.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.snr2.p24.t48.new.csv')
meta.s6.snr3.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.snr3.p24.t48.new.csv')
meta.s6.uneven.p24.t48.1.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.uneven.p24.t48.1.new.csv')
meta.s6.uneven.p24.t48.2.new <- read.csv('./sim_6_points_new/meta_rain/meta2d_s6.uneven.p24.t48.2.new.csv')



rain.s6.asy.p24.t48.f8.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.asy.p24.t48.f8.new.csv')
rain.s6.nonsta.p24.t48.f8.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.nonsta.p24.t48.f8.new.csv')
rain.s6.miss1.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/rain.s6.miss1.p24.t48.csv')
rain.s6.miss5.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/rain.s6.miss5.p24.t48.csv')
rain.s6.miss10.p24.t48 <- read.csv('./sim_6_points_new/meta_rain/rain.s6.miss10.p24.t48.csv')
rain.s6.snr0.5.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.snr0.5.p24.t48.new.csv')
rain.s6.snr1.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.snr1.p24.t48.new.csv')
rain.s6.snr2.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.snr2.p24.t48.new.csv')
rain.s6.snr3.p24.t48.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.snr3.p24.t48.new.csv')
rain.s6.uneven.p24.t48.1.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.uneven.p24.t48.1.new.csv')
rain.s6.uneven.p24.t48.2.new <- read.csv('./sim_6_points_new/meta_rain/rain.s6.uneven.p24.t48.2.new.csv')



write.table(s6.nonsta.p24.t48.f8.new[,-1], './sim_6_points_new/bc/data/s6.nonsta.p24.t48.f8.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.asy.p24.t48.f8.new[,-1], './sim_6_points_new/bc/data/s6.asy.p24.t48.f8.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.uneven.p24.t48.1.new[,-1], './sim_6_points_new/bc/data/s6.uneven.p24.t48.1.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.uneven.p24.t48.2.new[,-1], './sim_6_points_new/bc/data/s6.uneven.p24.t48.2.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.snr0.5.p24.t48.new[,-1], './sim_6_points_new/bc/data/s6.snr0.5.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.snr1.p24.t48.new[,-1], './sim_6_points_new/bc/data/s6.snr1.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.snr2.p24.t48.new[,-1], './sim_6_points_new/bc/data/s6.snr2.p24.t48.new.tsv',
            sep = '\t', col.names = NA)
write.table(s6.snr3.p24.t48.new[,-1], './sim_6_points_new/bc/data/s6.snr3.p24.t48.new.tsv',
            sep = '\t', col.names = NA)


bc.s6.asy.p24.t48.f8.new <- read.table('./sim_6_points_new/bc/results/bc_s6.asy.p24.t48.f8.new.tsv', header = T)
bc.s6.nonsta.p24.t48.f8.new <- read.table('./sim_6_points_new/bc/results/bc_s6.nonsta.p24.t48.f8.new.tsv', header = T)
bc.s6.snr0.5.p24.t48.new <- read.table('./sim_6_points_new/bc/results/bc_s6.snr0.5.p24.t48.new.tsv', header = T)
bc.s6.snr1.p24.t48.new <- read.table('./sim_6_points_new/bc/results/bc_s6.snr1.p24.t48.new.tsv', header = T)
bc.s6.snr2.p24.t48.new <- read.table('./sim_6_points_new/bc/results/bc_s6.snr2.p24.t48.new.tsv', header = T)
bc.s6.snr3.p24.t48.new <- read.table('./sim_6_points_new/bc/results/bc_s6.snr3.p24.t48.new.tsv', header = T)
bc.s6.uneven.p24.t48.1.new <- read.table('./sim_6_points_new/bc/results/bc_s6.uneven.p24.t48.1.new.tsv', header = T)
bc.s6.uneven.p24.t48.2.new <- read.table('./sim_6_points_new/bc/results/bc_s6.uneven.p24.t48.2.new.tsv', header = T)



write.table(s6.nonsta.p24.t48.f8.new[,-1], './sim_6_points_new/eJTK/s6.nonsta.p24.t48.f8.new.txt',
            sep = '\t', quote = F)
write.table(s6.asy.p24.t48.f8.new[,-1], './sim_6_points_new/eJTK/s6.asy.p24.t48.f8.new.txt',
            sep = '\t', quote = F)
write.table(s6.snr0.5.p24.t48.new[,-1], './sim_6_points_new/eJTK/s6.snr0.5.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s6.snr1.p24.t48.new[,-1], './sim_6_points_new/eJTK/s6.snr1.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s6.snr2.p24.t48.new[,-1], './sim_6_points_new/eJTK/s6.snr2.p24.t48.new.txt',
            sep = '\t', quote = F)
write.table(s6.snr3.p24.t48.new[,-1], './sim_6_points_new/eJTK/s6.snr3.p24.t48.new.txt',
            sep = '\t', quote = F)

eJTK.s6.asy.p24.t48.f8.new <- read.table('./sim_6_points_new/eJTK/s6.asy.p24.t48.f8.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s6.nonsta.p24.t48.f8.new <- read.table('./sim_6_points_new/eJTK/s6.nonsta.p24.t48.f8.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s6.snr0.5.p24.t48.new <- read.table('./sim_6_points_new/eJTK/s6.snr0.5.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s6.snr1.p24.t48.new <- read.table('./sim_6_points_new/eJTK/s6.snr1.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s6.snr2.p24.t48.new <- read.table('./sim_6_points_new/eJTK/s6.snr2.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)
eJTK.s6.snr3.p24.t48.new <- read.table('./sim_6_points_new/eJTK/s6.snr3.p24.t48.new_cos24_ph0-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt', header = T)


auc_6_table_new <- function(){
  ds <- AUC_gp1(meta.s6.nonsta.p24.t48.f8.new, rain.s6.nonsta.p24.t48.f8.new, bc.s6.nonsta.p24.t48.f8.new, eJTK.s6.nonsta.p24.t48.f8.new, 6000, 6000)
  ds <- rbind(ds, AUC_gp1(meta.s6.asy.p24.t48.f8.new, rain.s6.asy.p24.t48.f8.new, bc.s6.asy.p24.t48.f8.new, eJTK.s6.asy.p24.t48.f8.new, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp3(meta.s6.uneven.p24.t48.1.new, rain.s6.uneven.p24.t48.1.new, bc.s6.uneven.p24.t48.1.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp3(meta.s6.uneven.p24.t48.2.new, rain.s6.uneven.p24.t48.2.new, bc.s6.uneven.p24.t48.2.new, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp4(meta.s6.miss1.p24.t48, rain.s6.miss1.p24.t48, 6000, 6000))
  ds <- rbind(ds, AUC_gp4(meta.s6.miss5.p24.t48, rain.s6.miss5.p24.t48, 6000, 6000))
  ds <- rbind(ds, AUC_gp4(meta.s6.miss10.p24.t48, rain.s6.miss10.p24.t48, 6000, 6000))
  
  ds <- rbind(ds, AUC_gp1(meta.s6.snr0.5.p24.t48.new, rain.s6.snr0.5.p24.t48.new, bc.s6.snr0.5.p24.t48.new, eJTK.s6.snr0.5.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s6.snr1.p24.t48.new, rain.s6.snr1.p24.t48.new, bc.s6.snr1.p24.t48.new, eJTK.s6.snr1.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s6.snr2.p24.t48.new, rain.s6.snr2.p24.t48.new, bc.s6.snr2.p24.t48.new, eJTK.s6.snr2.p24.t48.new, 6000, 6000))
  ds <- rbind(ds, AUC_gp1(meta.s6.snr3.p24.t48.new, rain.s6.snr3.p24.t48.new, bc.s6.snr3.p24.t48.new, eJTK.s6.snr3.p24.t48.new, 6000, 6000))
  return(ds)
}

AUC_6_table_new <- auc_6_table_new()
AUC_6_table_new$AUC <- as.numeric(as.character(AUC_6_table_new$AUC))
write.csv(AUC_6_table_new, './sim_6_points_new/AUC_6_table_new.csv')




# AUC plot
library(tidyr)
auc.s6.group2.new <- matrix(nrow = 3, ncol = 7, dimnames = list(c('Stationary', 'Non-stationary','Asymmetric'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
auc.s6.group2.new[1,] <- AUC_6_table[c(1:7),3]
auc.s6.group2.new[2,] <- AUC_6_table_new[c(1:7),3]
auc.s6.group2.new[3,] <- AUC_6_table_new[c(8:14),3]
auc.s6.group2.new.ds <- gather(data.frame(auc.s6.group2.new, Category = row.names(auc.s6.group2.new),
                                      row.names = NULL), Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
auc.s6.group2.new.ds$Method <- factor(auc.s6.group2.new.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))
auc.s6.group2.new.ds$Category <- factor(auc.s6.group2.new.ds$Category, levels = c('Stationary','Non-stationary','Asymmetric'))



auc.s6.group3.new <- matrix(nrow = 2, ncol = 3, dimnames = list(c('1', '2'),c('LS/MC','RAIN','BC')))
auc.s6.group3.new[1,] <- AUC_6_table_new[c(15:17),3]
auc.s6.group3.new[2,] <- AUC_6_table_new[c(18:20),3]
auc.s6.group3.new.ds <- gather(data.frame(auc.s6.group3.new, Category = row.names(auc.s6.group3.new),
                                      row.names = NULL), Method, AUC,LS.MC,RAIN,BC)
auc.s6.group3.new.ds$Method <- factor(auc.s6.group3.new.ds$Method, levels = c('LS.MC','RAIN','BC'),
                                  labels = c('LS/MC','RAIN','BC'))


auc.s6.group4.new <- matrix(nrow = 3, ncol = 4, dimnames = list(c('1%', '5%', '10%'),c('LS','JTK','RAIN','MC')))
auc.s6.group4.new[1,] <- AUC_6_table_new[c(21:24),3]
auc.s6.group4.new[2,] <- AUC_6_table_new[c(25:28),3]
auc.s6.group4.new[3,] <- AUC_6_table_new[c(29:32),3]
auc.s6.group4.new.ds <- gather(data.frame(auc.s6.group4.new, Category = row.names(auc.s6.group4.new),
                                      row.names = NULL), Method, AUC,LS,JTK,RAIN,MC)
auc.s6.group4.new.ds$Method <- factor(auc.s6.group4.new.ds$Method, levels = c('LS','JTK','RAIN','MC'))
auc.s6.group4.new.ds$Category <- factor(auc.s6.group4.new.ds$Category, levels = c('1%', '5%', '10%'))


auc.s6.group6.new <- matrix(nrow = 4, ncol = 7, dimnames = list(c('SNR 0.5:1','SNR 1:1','SNR 2:1','SNR 3:1'),c('LS','ARS','JTK','RAIN','eJTK','MC','BC')))
auc.s6.group6.new[1,] <- AUC_6_table_new[c(33:39),3]
auc.s6.group6.new[2,] <- AUC_6_table_new[c(33:39)+7,3]
auc.s6.group6.new[3,] <- AUC_6_table_new[c(33:39)+14,3]
auc.s6.group6.new[4,] <- AUC_6_table_new[c(33:39)+21,3]
auc.s6.group6.new.ds <- gather(data.frame(auc.s6.group6.new, Category = row.names(auc.s6.group6.new),
                                      row.names = NULL), Method, AUC, LS, ARS, JTK, RAIN,eJTK,MC,BC)
auc.s6.group6.new.ds$Method <- factor(auc.s6.group6.new.ds$Method, levels = c('LS','ARS','JTK','RAIN','eJTK','MC','BC'))
auc.s6.group6.new.ds$Category <- factor(auc.s6.group6.new.ds$Category, levels = c('SNR 3:1', 'SNR 2:1', 'SNR 1:1', 'SNR 0.5:1'),
                                    labels = c('3:1', '2:1', '1:1', '0.5:1'))





# plot with legend
# group 1
library(ggplot2)
pdf(file = './figures_for_paper/AUC_simulation/6 points/6points_new_with_legend.pdf', height = 3, width = 4.5)
ggplot(data = auc.s4.group1.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.65,0.35)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') 
  

# group 2
ggplot(data = auc.s6.group2.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.30)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Waveforms', shape = 'Waveforms')

# group 3
ggplot(data = auc.s6.group3.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.33)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Uneven\ntime-points', shape = 'Uneven\ntime-points')

# group 4
ggplot(data = auc.s6.group4.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.33)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'Missing values', shape = 'Missing values')

# group 5
ggplot(data = auc.s4.group5.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  #theme(legend.position = c(0.7,0.35)) +
  #theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Replicates', shape = 'Replicates')

# group 6
ggplot(data = auc.s6.group6.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  # theme(legend.position = c(0.7,0.29)) +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  # theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'SNR', shape = 'SNR') + guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2))
dev.off()




# plot without legend
pdf(file = './figures_for_paper/AUC_simulation/6 points/6points_new_without_legend.pdf', height = 3, width = 3.5)
ggplot(data = auc.s4.group1.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Sampling patterns', shape = 'Sampling patterns') +
  theme(panel.border = element_blank())


# group 2
ggplot(data = auc.s6.group2.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Waveforms', shape = 'Waveforms') +
  theme(panel.border = element_blank())

# group 3
ggplot(data = auc.s6.group3.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Uneven\ntime-points', shape = 'Uneven\ntime-points') +
  theme(panel.border = element_blank())

# group 4
ggplot(data = auc.s6.group4.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) +
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'Missing values', shape = 'Missing values') +
  theme(panel.border = element_blank())

# group 5
ggplot(data = auc.s4.group5.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  #theme(axis.title = element_text(size = rel(1.1)), axis.text = element_text(size = rel(1.1))) + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = 'Replicates', shape = 'Replicates') +
  theme(panel.border = element_blank())

# group 6
ggplot(data = auc.s6.group6.new.ds, aes(x = Method, y = AUC, color = Category, shape = Category, group = Category)) + geom_point() + 
  geom_line() + ylim(0.3,1) + theme_bw() + theme_bw(base_size=12) + 
  theme(legend.position = 'none') +
  # theme(axis.title = element_text(size = rel(1.3)), axis.text = element_text(size = rel(1.1))) + 
  # theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(color = 'SNR', shape = 'SNR') + guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2)) + 
  theme(panel.border = element_blank())
  
dev.off()




