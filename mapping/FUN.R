B.plot <- function(B, r = range(B), ...) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)

  
  plot(NA, ylim = c(1, n), xlim = r, axes=FALSE, ...)
  axis(side =1)
  axis(side =2, at = 1:n, labels = colnames(B), las = 2)
  
  for(i in 1:n) {
  	lines(B[c(1, 5),i], c(i,i))
  	lines(B[c(2, 4),i], c(i,i), lwd = 3)
  	points(B[3,i], i, pch = 20, cex = 1, col = 2)
  }

}
B.plotAdd <- function(B, offset = .5, col = 4) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)
  r <- range(B)
    
  for(i in 1:n) {
  	lines(B[c(1, 5),i], c(i,i) + offset, col = col)
  	lines(B[c(2, 4),i], c(i,i) + offset, lwd = 1.5, col = col)
  	points(B[3,i], i+offset, pch = 20, cex = 1, col = col)
  }

}

B.plotV <- function(B, ...) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)
  r <- range(B)
  
  plot(NA, xlim = c(1, n), ylim = r, axes=FALSE, ...)
  axis(side =2)
  axis(side =1, at = 1:n, labels = colnames(B), las = 2)
  
  for(i in 1:n) {
  	lines(c(i,i), B[c(1, 5),i])
  	lines(c(i,i), B[c(2, 4),i], lwd = 3)
  	points(i, B[3,i], pch = 20, cex = 1, col = 2)
  }

}

B.plotXV <- function(B, x, ...) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)
  r <- range(B)
  
  plot(NA, xlim = c(0, max(x)), ylim = r, ...)
  for(i in 1:n) {
  	lines(x[c(i,i)], B[c(1, 5),i])
  	lines(x[c(i,i)], B[c(2, 4),i], lwd = 3)
  	points(x[i], B[3,i], pch = 20, cex = 1, col = 2)
  }

}

B.plotVAdd <- function(B, offset = .5) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)
  r <- range(B)
    
 for(i in 1:n) {
  	lines(c(i,i) + offset, B[c(1, 5),i], col = 4)
  	lines(c(i,i) + offset, B[c(2, 4),i], lwd = 1.5, col = 4)
  	points(i + offset, B[3,i], pch = 20, cex = 1, col = 4)
  }

}

logsum <- function(x) {
  M <- max(x)
  M + log(sum(exp(x-M)))
}
colVars <- function(x) apply(x, 2, var)
waic <- function(log_lik) {
  
  ## Handle the dimensions
  if(length(dim(log_lik))==1) {
  	dim(log_lik) <- c(length(log_lik), 1)
  } else {
  	c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  }
  
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  
  ### Compute the waic
  lpd <- apply(log_lik, 2, logsum) - log(S)  # log(mean(p(y_i|theta^s)))
  p_waic <- colVars(log_lik)                 # var(log(p(y_i|theta^s)))
  p_waic_2 <- 2 * (lpd - colMeans(log_lik))
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  
  ### Compute the LOO
  loo_weights_raw <- 1/exp(log_lik - max(log_lik))
  loo_weights_normalized <- loo_weights_raw/matrix(colMeans(loo_weights_raw),nrow = S, ncol=n,byrow=TRUE)
  l_loo_weights_regularized <- log(pmin(loo_weights_normalized, sqrt(S)))
  elpd_loo <- apply(log_lik + l_loo_weights_regularized, 2, logsum) - apply(l_loo_weights_regularized, 2, logsum)
  p_loo <- lpd - elpd_loo
  pointwise <- cbind('waic' = waic, 'lpd' = lpd, 'p_waic' = p_waic, 'p_w2' = p_waic_2, 'elpd_waic' = elpd_waic, 
    'p_loo' = p_loo, 'elpd_loo' = elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n * colVars(pointwise))
  ### 
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
    p_waic=total["p_waic"], p_waic2 = total["p_w2"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
    pointwise=pointwise, total=total, se=se))
}



qProb <- function(x, probs = c(.025, .158, .5, 1-.158, .975)) {
  c(quantile(x, probs = probs), 'M' = mean(x), 'V' = var(x))
}

## XL contains list of (matirx of values, Multiplier (for diet), and an extender.. )
VarExp <- function(X_L, Den, w = 1) {
  if(!is.matrix(Den)) {
  	DV <- var(Den*w)
  } else {
  	DV <- apply(Den, 1, function(x) var(x * w))
  }
  
  Out <- sapply(X_L, function(x) {
  	qProb(apply(x[[1]], 1, function(xi) var((xi[x[[3]]]) * x[[2]] ))/DV)
  })
  
  
  Out
}


VarExpQTL_main <- function(X_L, Den, w = 1) {
  if(!is.matrix(Den)) {
    DV <- var(Den*w)
  } else {
    DV <- apply(Den, 1, function(x) var(x * w))
  }
  
  ## Need to substract variances so don't summarize with qProb yet
  Out <- sapply(X_L, function(x) {
    apply(x[[1]], 1, function(xi) var((xi[x[[3]]]) * x[[2]] ))/DV
  }, simplify = FALSE)
  
  A <- qProb(Out[['A']])
  DE <- qProb(Out[['D']] - Out[['A']])
  D <- qProb(Out[['D']])
  FE <- qProb(Out[['QTL']] - Out[['D']])
  QTL <- qProb(Out[['QTL']])
  cbind(A, DE, D, FE, QTL)
}


VarExpQTL_diet <- function(X_L, Den, w = 1) {
  if(!is.matrix(Den)) {
    DV <- var(Den*w)
  } else {
    DV <- apply(Den, 1, function(x) var(x * w))
  }
  
  ## Need to substract variances so don't summarize with qProb yet
  Out <- sapply(X_L, function(x) {
    apply(x[[1]], 1, function(xi) var((xi[x[[3]]]) * x[[2]] ))/DV
  }, simplify = FALSE)
  
  AXD <- qProb(Out[['AXD']])
  DEXD <- qProb(Out[['DXD']] - Out[['AXD']])
  DXD <- qProb(Out[['DXD']])
  FEXD <- qProb(Out[['QTLXD']] - Out[['DXD']])
  QTLXD <- DXD + FEXD
  cbind(AXD, DEXD, DXD, FEXD, QTLXD)
}




R2Lam <- function(m, th) {
  if(!is.matrix(th)) {
  	dR2 <- var(th)
  	e <- sweep(-m, 2, th, FUN = '+')
  } else {
  	dR2 <- mean(apply(th, 1, var))
  	e <- th - m
  }
  nR2 <- mean(apply(e, 1, var))
  nL <- var(colMeans(e))
  
  c(1-nR2/dR2, 1-nL/nR2)
}

fLCr <- function(x,M,P) x[,M] + x[,P]
fLCrG <- function(x,M,P) x[,M] - x[,P]
## This one is special because there are 2 Lines
## QTL comes w/ an X in the first and the way to get to the next in the others..
## I.e. X[[1]] * B[[1]], X[[1]] * X[[2]] * B[[2]], X [[1]] * 
## L is the gorup of L's and w/ 2 groups.. 
VCrExp <- function(QTL = NA, L, Cr) {
  DV <- apply(Cr, 1, var)  
  qtpP <- numeric(0)
  if(!is.na(QTL[1])) {
  	qtpP <- sapply(QTL, function(x) qProb(apply(x[[1]], 1, var)/DV))
  }  

  L <- qProb(apply(L[[1]][,L[[2]]] + L[[1]][,L[[3]]], 1, var)/DV)  
  cbind(qtpP, L)  
}

VCrExpG <- function(QTL = NA, L, G, Cr) {
  qtpP <- numeric(0)
  if(!is.na(QTL[1])) {
  	
  }  
  DV <- apply(Cr, 1, var)  
  L <- qProb(apply(L[[1]][,L[[2]]] + L[[1]][,L[[3]]], 1, var)/DV)  
  G <- qProb(apply(G[[1]][,G[[2]]] - G[[1]][,G[[3]]], 1, var)/DV)    
  cbind(qtpP, L, G)  
}


## This function is about running the code for A/D/F/N all together..
## L will be a lsit of things where first component is [A,D] matrix and there are others
## used for each part..
## Cr is in there
## f handles the rest.. 
AD.Split <- function(L, f) {
  LA <- LD <- LN <- LF <- L
  for(i in 1:length(L)) {
  	LA[[i]][[1]] <- LA[[i]][[1]][,1,] 
  	LD[[i]][[1]] <- LD[[i]][[1]][,2,] 
  	LN[[i]][[1]] <- LN[[i]][[1]][,1,] - .5 * LN[[i]][[1]][,2,]
  	LF[[i]][[1]] <- LF[[i]][[1]][,1,] + .5 * LF[[i]][[1]][,2,]
  }
  
  list('A' = f(LA), 'D' = f(LD), 'N' = f(LN), 'F' = f(LF))
}


B.plot2 <- function(B, ...) { ## assume a matrix of quantiles.. colnames are names, rownaems are quants..
  n <- ncol(B)
  p <- nrow(B)
  r <- range(B)
  
  o <- order(B[3,])
  
  plot(NA, xlim = c(1, n), ylim = r, axes=FALSE, ...)
  axis(side =2)
  axis(side =1, at = 1:n, labels = colnames(B))

  lines(B[3,o])
  lines(B[1,o], col = 2, lty = 3)
  lines(B[5,o], col = 2, lty = 3)
  lines(B[2,o], col = 3, lty = 2)
  lines(B[4,o], col = 3, lty = 2)
 
}