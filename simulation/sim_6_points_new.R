#! /usr/bin/Rscript

## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)



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



miss_ing2 <- function(rate, dataset, miss_col = c(4,9,12)){
  nrow <- dim(dataset)[1]
  miss_row <- sample(x = c(1:nrow), size = nrow*rate/100)
  dataset[miss_row, miss_col] <- NA
  return(dataset)
}


# measure_sequence <- function(colnames, n){
#   numbers <- as.numeric(gsub("CT", "", colnames[-1]))
#   time_seq <- seq(1,48,1)
#   rslt <- rep(0,48)
#   for (i in 1:n)
#     rslt[numbers[i]+1] <- 1
#   return(list(results = rslt, numbers = numbers))
# }

measure_sequence <- function(colnames, n, freq){
  numbers <- as.numeric(gsub("CT", "", colnames[-1]))
  rslt <- rep(0,n)
  for (i in 1:n)
    rslt[numbers[i]/freq+1] <- 1
  return(list(results = rslt, numbers = numbers))
}


runit <- function(object = s3.nonsta.p24.t48.f4.new, time = seq(0,44,4), freq = 4, time2 = NULL){
  library(MetaCycle)
  library(rain)
  
  s1 <- paste(paste('/',deparse(substitute(object)),sep = ''), '.csv', sep = '')
  s2 <- paste(paste('/rain.',deparse(substitute(object)),sep = ''), '.csv', sep = '')
  
  write.csv(object[,-1], paste(args$indir, s1, sep = ''))
  meta2d(infile = paste(args$indir, s1, sep = ''),
         outputFile = T,  filestyle = 'csv',
         timepoints = time, outIntegration = 'onlyIntegration', maxper = 28,
         minper = 20, outdir = args$outdir, parallelize = T, nCores = 4)
  
  if (is.null(time2)) {
    rain.s3.nonsta.p24.t48.f4.new <- rain(t(object[,-1]), deltat=freq, period=24, period.delta = 4,
                                          peak.border = c(0.3, 0.7), nr.series = 1, adjp.method = "ABH")
  }
  else {
    rain.s3.nonsta.p24.t48.f4.new <- rain(t(object[,-1]), deltat=freq, period=24, period.delta = 4,
                                          peak.border = c(0.3, 0.7), nr.series = 1, adjp.method = "ABH",
                                          measure.sequence = time2)
  }
  write.csv(rain.s3.nonsta.p24.t48.f4.new, paste(args$outdir, s2, sep = ''))
}







# simulation & analysis
# group1
set.seed(101)
s6.p24.t24.f2.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                          freq = 2, type=c('cos','cos2','cos_peak','flat'))
s6.p24.t48.f4.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 47,
                          freq = 4, type=c('cos','cos2','cos_peak','flat'))
s6.p24.t48.f8.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 47,
                          freq = 8, type=c('cos','cos2','cos_peak','flat'))
s6.p24.t24.f4.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                          freq = 4, type=c('cos','cos2','cos_peak','flat'))
runit(s6.p24.t24.f2.new, time = seq(0,22,2), freq = 2)
runit(s6.p24.t48.f4.new, time = seq(0,44,4), freq = 4)
runit(s6.p24.t48.f8.new, time = seq(0,40,8), freq = 8)
runit(s6.p24.t24.f4.new, time = seq(0,20,4), freq = 4)




# group2
# stastionary is the same as the s3.p24.t48.f4.new
set.seed(102)
s6.nonsta.p24.t48.f8.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23, freq = 4,
                                 type=c('cos_damp','trend_exp', 'trend_linear','linear'))
s6.asy.p24.t48.f8.new <- comb(n = 2, amount = c(6000,6000), time_window = 23, freq = 4,
                             type=c('saw_tooth','flat'))
runit(s6.nonsta.p24.t48.f8.new, time = seq(0,20,4), freq = 4)
runit(s6.asy.p24.t48.f8.new, time = seq(0,20,4), freq = 4)


# new.group3
set.seed(103)
s6.uneven.p24.t48.1.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                                freq = 4, type=c('cos','cos2','cos_peak','flat'))[,-4]
s6.uneven.p24.t48.2.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                                freq = 4, type=c('cos','cos2','cos_peak','flat'))[,-c(3,5)]
runit(object = s6.uneven.p24.t48.1.new,
      time = measure_sequence(colnames = colnames(s6.uneven.p24.t48.1.new), n = 6, freq = 4)$numbers,
      time2 = measure_sequence(colnames = colnames(s6.uneven.p24.t48.1.new), n = 6, freq = 4)$results)
runit(object = s6.uneven.p24.t48.2.new,
      time = measure_sequence(colnames = colnames(s6.uneven.p24.t48.2.new), n = 6, freq = 4)$numbers,
      time2 = measure_sequence(colnames = colnames(s6.uneven.p24.t48.2.new), n = 6, freq = 4)$results)




# new.group4
set.seed(104)
s6.miss1.p24.t48 <- miss_ing2(dataset = comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                                   freq = 4, type=c('cos','cos2','cos_peak','flat')),
                              rate = 1, miss_col = c(3,5,6))
s6.miss5.p24.t48 <- miss_ing2(dataset = comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                                   freq = 4, type=c('cos','cos2','cos_peak','flat')),
                              rate = 5, miss_col = c(3,5,6))
s6.miss10.p24.t48 <- miss_ing2(dataset = comb(n = 4, amount = c(rep(2000,3),6000), time_window = 23,
                                    freq = 4, type=c('cos','cos2','cos_peak','flat')),
                               rate = 10, miss_col = c(3,5,6))

runit(s6.miss1.p24.t48, time = seq(0,20,4), freq = 4)
runit(s6.miss5.p24.t48, time = seq(0,20,4), freq = 4)
runit(s6.miss10.p24.t48, time = seq(0,20,4), freq = 4)


# new.group5
set.seed(105)
s6.r2.p24.t48.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 95,
                          freq = 8, type=c('cos','cos2','cos_peak','flat'))
colnames(s6.r2.p24.t48.new)[-1] <- paste(rep(paste('CT', seq(0,40,8), sep = ''),2), paste('_',c(sapply(seq(1,2,1),function(x){rep(x,6)})), sep=''), sep = '')

s6.r2.p24.t24.new <- comb(n = 4, amount = c(rep(2000,3),6000), time_window = 47,
                          freq = 8, type=c('cos','cos2','cos_peak','flat'))
colnames(s6.r2.p24.t24.new)[-1] <- paste(rep(paste('CT', seq(0,16,8), sep = ''),2), paste('_',c(sapply(seq(1,2,1),function(x){rep(x,3)})), sep=''), sep = '')

runit(s6.r2.p24.t48.new, time = rep(seq(0,40,8),2), freq = 8)
runit(s6.r2.p24.t24.new, time = rep(seq(0,16,8),2), freq = 8)



# new.group6
set.seed(106)
s6.snr0.5.p24.t48.new <- comb(n = 2, amount = c(6000,6000), time_window = 23, freq = 4,
                              type=c('SNR_cos','flat'), SNR = 0.5)
s6.snr1.p24.t48.new <- comb(n = 2, amount = c(6000,6000), time_window = 23, freq = 4,
                            type=c('SNR_cos','flat'), SNR = 1)
s6.snr2.p24.t48.new <- comb(n = 2, amount = c(6000,6000), time_window = 23, freq = 4,
                            type=c('SNR_cos','flat'), SNR = 2)
s6.snr3.p24.t48.new <- comb(n = 2, amount = c(6000,6000), time_window = 23, freq = 4,
                            type=c('SNR_cos','flat'), SNR = 3)

runit(s6.snr0.5.p24.t48.new, time = seq(0,20,4), freq = 4)
runit(s6.snr1.p24.t48.new, time = seq(0,20,4), freq = 4)
runit(s6.snr2.p24.t48.new, time = seq(0,20,4), freq = 4)
runit(s6.snr3.p24.t48.new, time = seq(0,20,4), freq = 4)




