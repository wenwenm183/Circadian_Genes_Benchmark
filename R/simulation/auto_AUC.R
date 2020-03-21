# There are two functions that are both suitable for all of the methods.
# The first V_auc4 is designed for 4-shape dataset, and the second is designed
# for 2-shape. For every time, please only input ONE method. The input ds is the 
# result of each method, don't need to do any change on it.

auto_auc <- function(ds, method, n.p, n.n){
  library(qvalue)
  library(pROC)
  
  cl <- switch (method,
                'ARSER' = 'ARS_BH.Q', 'JTK' = 'JTK_BH.Q', 'LS' = 'LS_BH.Q', 
                'meta2d' = 'meta2d_BH.Q', 'rain' = 'pVal', 'bc' = 'Q_VALUE', 
                'eJTK' = 'GammaBH')
  if (method != 'rain') {
    AUC <- suppressMessages(roc(c(rep(1,n.p),rep(0,n.n)), ds[,cl]))$auc
  }
  else {
    Qvalue <- qvalue(ds[,cl])$qvalue
    AUC <- suppressMessages(roc(c(rep(1,n.p),rep(0,n.n)), Qvalue))$auc
  }
  return(list(method = method, AUC = AUC))
}



AUC_gp1 <- function(meta, rain, bc, eJTK, n.p, n.n){
    ds <- as.data.frame(do.call("cbind", auto_auc(meta, method = 'LS',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'ARSER',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'JTK',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(rain, method = 'rain',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(eJTK, method = 'eJTK',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'meta2d',n.p, n.n)))
    ds <- rbind(ds, do.call("cbind", auto_auc(bc, method = 'bc',n.p, n.n)))
    ds$group <- rep(unlist(strsplit(deparse(substitute(meta)), 'meta.'))[2], 7)
    ds <- ds[,c(3,1,2)]
    return(ds)
}




AUC_gp3 <- function(meta, rain, bc, n.p, n.n){
  ds <- as.data.frame(do.call("cbind", auto_auc(meta, method = 'LS',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(rain, method = 'rain',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(bc, method = 'bc',n.p, n.n)))
  ds$group <- rep(unlist(strsplit(deparse(substitute(meta)), 'meta.'))[2], 3)
  ds <- ds[,c(3,1,2)]
  return(ds)
}


AUC_gp4 <- function(meta, rain, n.p, n.n){
  ds <- as.data.frame(do.call("cbind", auto_auc(meta, method = 'LS',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'JTK',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(rain, method = 'rain',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'meta2d',n.p, n.n)))
  ds$group <- rep(unlist(strsplit(deparse(substitute(meta)), 'meta.'))[2], 4)
  ds <- ds[,c(3,1,2)]
  return(ds)
}


AUC_gp5 <- function(meta, rain, bc, eJTK, n.p, n.n){
  ds <- as.data.frame(do.call("cbind", auto_auc(meta, method = 'LS',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'JTK',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(rain, method = 'rain',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(eJTK, method = 'eJTK',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(meta, method = 'meta2d',n.p, n.n)))
  ds <- rbind(ds, do.call("cbind", auto_auc(bc, method = 'bc',n.p, n.n)))
  ds$group <- rep(unlist(strsplit(deparse(substitute(meta)), 'meta.'))[2], 6)
  ds <- ds[,c(3,1,2)]
  return(ds)
}



