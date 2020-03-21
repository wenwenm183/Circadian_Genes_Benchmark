# this function generates one piece of time series (i.e. one sample). It is the
# input of simulation()
one_row <- function(per = 24, time_window, freq, type, SNR){
  pha <- runif(1, 0, per)
  noise <- rnorm((time_window+1)/freq, 0, 1)
  amp <- runif(1, 1, 6)
  t <- seq(0, time_window, freq)
  # xt <- amp*cos(2*pi/per*(t - pha)) + rnorm(time_window+1, 0, 1)
  per2 <- 1/3*per
  pha1 <- pha+0.215*per/pi
  pha2 <- (pha1+0.25*per2)%%per
  slope <- runif(1, -0.05, 0)
  xt <- switch (type,
    'cos' = amp*cos(2*pi/per*(t - pha)) + noise,
    'cos2' = amp/1.39*(cos(2*pi/per*(t - pha1) ) + 
                         0.5*cos( 2*pi/per2*(t - pha2))) + noise,
    'cos_peak' = amp*(-1 + 2*abs(cos(pi/per*(t - pha)))**10) + noise,
    'flat' = rep(0,(time_window+1)/freq) + noise,
    'linear' = slope*t + noise,
    'cos_damp' = amp*cos(2*pi/per*(t - pha))*exp(-0.01*t) + noise,
    'saw_tooth' = -2*(amp/pi)*atan(1/tan((pi/per)*(t - pha))) + noise,
    'trend_exp' = 5*exp(-0.01*t) + amp*cos(2*pi/per*(t - pha)) + noise,
    'trend_linear' = amp*cos(2*pi/per*(t - pha)) + slope*t + noise,
    'SNR_cos' = sqrt(SNR*2)*cos(2*pi/per*(t - pha)) + noise
  )
  return(xt)
}

# test one_row()
# one_row(per = 24, time_window = 23, freq = 2, type = 'linear')
# one_row(per = 24, time_window = 47, freq = 4, type = 'flat')
#one_row(per = 24, time_window = 47, freq = 4, type = 'cos')
#one_row(per = 24, time_window = 47, freq = 4, type = 'cos_peak')
#one_row(per = 24, time_window = 47, freq = 4, type = 'cos_damp')
#plot(seq(0,47,2), one_row(per = 24, time_window = 47, freq = 2, type = 'saw_tooth'))
#plot(seq(0,47,2), one_row(per = 24, time_window = 47, freq = 2, type = 'trend_exp'))
#plot(seq(0,47,2), one_row(per = 24, time_window = 47, freq = 2, type = 'trend_linear'))




# This function generates a dataframe of simulation. It is the input of comb()
simulation <- function(per = 24, time_window, freq, amount, type, SNR){
  smlt <- t(apply(matrix(nrow = amount, ncol = (time_window+1)/freq), 1, 
                  function(x){one_row(per,time_window, freq, type, SNR)}))
  smlt <- as.data.frame(smlt)
  colnames(smlt) <- paste('CT',seq(0,time_window, freq), sep = '')
  return(smlt)
}



# This function is able to generate simulations of 5 types. 
# n means the number of types of curves you want.
# amount means the amount of simulations of each type of curve.
# type specifies the name of the curve.
comb <- function(per = 24, time_window, freq, n, amount, type, SNR){
  dataset <- data.frame(type = rep(type[1], amount[1]), 
                        simulation(per, time_window, freq, amount = amount[1], 
                                   type = type[1], SNR))
  for (i in 2:n){
    dataset <- rbind(dataset, data.frame(type = rep(type[i], amount[i]), 
                                         simulation(per, time_window, 
                                                    freq, amount = amount[i], 
                                                    type = type[i], SNR)))
  }
  return(dataset)
}


miss_ing <- function(nrow, ncol, rate, dataset){
  index <- sample(1:(nrow*ncol), rate*nrow*ncol)
  row_index <- floor(index/ncol)+1
  col_index <- index%%ncol + 1
  missing.mtrx <- matrix(c(row_index,col_index), ncol = 2)
  for (i in 1:(rate*nrow*ncol)){
    dataset[missing.mtrx[i,][1], missing.mtrx[i,][2]] <- NA
  }
  return(dataset)
}



miss_ing2 <- function(rate, dataset){
  nrow <- dim(dataset)[1]
  miss_row <- sample(x = c(1:nrow), size = nrow*rate/100)
  miss_col <- c(4,9,12)
  dataset[miss_row, miss_col] <- NA
  return(dataset)
}









# Group1 for different shapes and time windows 
set.seed(101)
s.p24.t48.f1 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                     freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))
s.p24.t48.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                     freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
s.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                     freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
s.p24.t24.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                     freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
s.p24.t24.f1 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                     freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))

s.p24.t48.f1.b <- cbind(type = s.p24.t48.f1[, 1], as.data.frame(t(apply(s.p24.t48.f1[,-1], 1, function(x){x+runif(1,5,10)}))))
save(s.p24.t48.f1.b, file ='./simulation/sampling_patterns/s.p24.t48.f1.b.rda')

s.p24.t24.f2.b <- cbind(type = s.p24.t24.f2[, 1], as.data.frame(t(apply(s.p24.t24.f2[,-1], 1, function(x){x+runif(1,10,20)}))))
save(s.p24.t24.f2.b, file ='./simulation/sampling_patterns/s.p24.t24.f2.b.rda')

save(s.p24.t48.f1, file ='./simulation/sampling_patterns/s.p24.t48.f1.rda')
save(s.p24.t48.f2, file ='./simulation/sampling_patterns/s.p24.t48.f2.rda')
save(s.p24.t48.f4, file ='./simulation/sampling_patterns/s.p24.t48.f4.rda')
save(s.p24.t24.f2, file ='./simulation/sampling_patterns/s.p24.t24.f2.rda')
save(s.p24.t24.f1, file ='./simulation/sampling_patterns/s.p24.t24.f1.rda')

which(s.p24.t48.f1 < 0)


# Group2 for different shapes
set.seed(102)
s.shapes.p24.t48.f4 <- comb(n = 9, amount = c(rep(200,9)), time_window = 47, freq = 4,
                            type=c('cos','cos2','cos_peak','cos_damp','trend_exp',
                                   'trend_linear','saw_tooth','flat','linear'))
s.shapes.p24.t23.f2 <- comb(n = 9, amount = c(rep(200,9)), time_window = 23, freq = 2,
                            type=c('cos','cos2','cos_peak','cos_damp','trend_exp',
                                   'trend_linear','saw_tooth','flat','linear'))
save(s.shapes.p24.t48.f4, file ='./simulation/sampling_patterns/s.shapes.p24.t48.f4.rda')
save(s.shapes.p24.t23.f2, file ='./simulation/sampling_patterns/s.shapes.p24.t23.f2.rda')



# Group3 Uneven sampling
set.seed(103)
s.uneven.p24.t48.24 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                         freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))[,c(1,sort(sample(2:49,24)))]
s.uneven.p24.t48.30 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                         freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))[,c(1,sort(sample(2:49,30)))]
s.uneven.p24.t48.18 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                            freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))[,c(1,sort(sample(2:49,18)))]
s.uneven.p24.t48.12 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                            freq = 1, type=c('cos','cos2','cos_peak','flat','linear'))[,c(1,sort(sample(2:49,12)))]

save(s.uneven.p24.t48.24, file ='./simulation/sampling_patterns/s.uneven.p24.t48.24.rda')
save(s.uneven.p24.t48.30, file ='./simulation/sampling_patterns/s.uneven.p24.t48.30.rda')
save(s.uneven.p24.t48.18, file ='./simulation/sampling_patterns/s.uneven.p24.t48.18.rda')
save(s.uneven.p24.t48.12, file ='./simulation/sampling_patterns/s.uneven.p24.t48.12.rda')


# Group4 Missing Value (1%, 5%, 10%)
set.seed(104)
s.missing1.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                              freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing1.p24.t48.f4[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.01, dataset = s.missing1.p24.t48.f4[,-1])

s.missing5.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                              freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing5.p24.t48.f4[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.05, dataset = s.missing5.p24.t48.f4[,-1])

s.missing10.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                              freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing10.p24.t48.f4[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.1, dataset = s.missing10.p24.t48.f4[,-1])

s.missing1.p24.t24.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                              freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing1.p24.t24.f2[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.01, dataset = s.missing1.p24.t24.f2[,-1])

s.missing5.p24.t24.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                              freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing5.p24.t24.f2[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.05, dataset = s.missing5.p24.t24.f2[,-1])

s.missing10.p24.t24.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                               freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
s.missing10.p24.t24.f2[,-1] <- miss_ing(nrow = 1000, ncol = 12, rate = 0.1, dataset = s.missing10.p24.t24.f2[,-1])

save(s.missing1.p24.t48.f4, file ='./simulation/sampling_patterns/s.missing1.p24.t48.f4.rda')
save(s.missing5.p24.t48.f4, file ='./simulation/sampling_patterns/s.missing5.p24.t48.f4.rda')
save(s.missing10.p24.t48.f4, file ='./simulation/sampling_patterns/s.missing10.p24.t48.f4.rda')
save(s.missing1.p24.t24.f2, file ='./simulation/sampling_patterns/s.missing1.p24.t24.f2.rda')
save(s.missing5.p24.t24.f2, file ='./simulation/sampling_patterns/s.missing5.p24.t24.f2.rda')
save(s.missing10.p24.t24.f2, file ='./simulation/sampling_patterns/s.missing10.p24.t24.f2.rda')


# Group5 Replicates (1 vs. 2)
set.seed(105)
s.r2.p24.t48 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 95,
                     freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
plot(seq(1,24,1), s.r2.p24.t48[1,-1])

s.r2.p24.t24 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                     freq = 2, type=c('cos','cos2','cos_peak','flat','linear'))
colnames(s.r2.p24.t48)[-1] <- paste(rep(paste('CT', seq(0,44,4), sep = ''),2), paste('_',c(sapply(seq(1,2,1),function(x){rep(x,12)})), sep=''), sep = '')
colnames(s.r2.p24.t24)[-1] <- paste(rep(paste('CT', seq(0,22,2), sep = ''),2), paste('_',c(sapply(seq(1,2,1),function(x){rep(x,12)})), sep=''), sep = '')

save(s.r2.p24.t48, file ='./simulation/sampling_patterns/s.r2.p24.t48.rda')
save(s.r2.p24.t24, file ='./simulation/sampling_patterns/s.r2.p24.t24.rda')



# new group for detrending
set.seed(1021)
s.detrend.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 47,
                             freq = 4, type=c('cos','trend_exp','trend_linear','flat','linear'))

s.detrend.p24.t24.f2 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 23,
                             freq = 2, type=c('cos','trend_exp','trend_linear','flat','linear'))



# dataset for replicates percentage missing values
set.seed(1051)
s.r3.missing10.p24.t48.f4 <- comb(n = 5, amount = c(200,200,200,200,200), time_window = 48*3-1,
                                  freq = 4, type=c('cos','cos2','cos_peak','flat','linear'))
s.r3.missing10.p24.t48.f4[,-1] <- miss_ing(nrow = 1000, ncol = 36, rate = 0.1, dataset = s.r3.missing10.p24.t48.f4[,-1])
colnames(s.r3.missing10.p24.t48.f4)[-1] <- paste(rep(paste('CT', seq(0,44,4), sep = ''),3), paste('_',c(sapply(seq(1,3,1),function(x){rep(x,12)})), sep=''), sep = '')

aa <- s.r3.missing10.p24.t48.f4[c(1:10),-1]
aa[1,1] <- NA

aaaa = data.frame(id = 1, CT0_1 = NA, CT1_1 = 1, CT2_1 = 1, 
                 CT0_2 = 1, CT1_2 = NA, CT2_2 = 1,
                 CT0_3 = 1, CT1_3 = 1, CT2_3 = NA)
write.table('./')









# Final simulation 1000 per shape
## Group5 Replicates
set.seed(205)
s.

