MM_BlockDiag <- function(y, X, Zt, nZ) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)
  Lam <- Diagonal(nrow(Zt), 2)
  mapLam <- rep(1:length(nZ), nZ)
  L0 <- Cholesky(tcrossprod(Zt) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
	
  dev <- function(theta) {
    Lam <- Diagonal(nrow(Zt), sqrt(theta)[mapLam])  ## specific to this model..  
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, -crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), -Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	# print(r2)
## Multiply by 2 is part of this.. 
	as.numeric(2 * determinant(L)$modulus + n * (1 + log(2 * pi * r2/n)))
  }
	
  cdev <- cmpfun(dev)
  cdev
}
DiACr <- cmpfun(MM_BlockDiag)


### Grad w/ s2 profiled out.. This is odne to help w/ stability for diagonal MM's but
### it does tend to be slower.. It is nice that it isn't necessary.. 
MM_BlockDiagGG <- function(y, X, Zt, nZ) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)
  Lam <- Diagonal(nrow(Zt), 2)
  mapLam <- rep(1:length(nZ), nZ)
  L0 <- Cholesky(tcrossprod(Zt) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
  w <- rep(1:length(nZ), nZ)
  dev <- function(theta) {
    Lam <- Diagonal(nrow(Zt), exp(theta)[mapLam])  ## specific to this model..  
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	l1 <- as.numeric(2 * determinant(L)$modulus)
    # # print(c(r2, l1, n * (1 + log(2 * pi * r2/n))))

    ## Multiply by 2 is part of this.. 
	c(l1 + n * (1 + log(2 * pi * r2/n)), r2, l1, n * (1 + log(2 * pi * r2/n)))
  }
  
  grad <- function(theta) {
    Lam <- Diagonal(nrow(Zt), exp(theta)[mapLam])  ## specific to this model..  
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- backsolve(RX, Xty - crossprod(RXZ, cu), transpose = TRUE)
  	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
  	B.hat <- backsolve(RX, cB)
  	u.hat <- solve(L, solve(L, cu - RXZ %*% B.hat, system = 'Lt'), 'Pt')
  	   	 
  	DL <- solve(L, solve(L, Diagonal(nrow(Zt), 1), system = 'P'), 'L')
    th1 <- colSums(DL * DL)
    th1 <- 2 * (nZ - as.numeric(tapply(th1, w, sum)))
    th2 <- 2*n/r2 * as.numeric(tapply(u.hat * u.hat, w, sum))
    th1 <<- th1
    th2 <<- th2

	th1 -  th2
  }

  cholT <- function(theta) {
  	Lam <- Diagonal(nrow(Zt), exp(theta)[mapLam])
  	ZtP <- crossprod(Lam, Zt)
  	L <- update(L0, ZtP, mult = 1)
  	LV <- Cholesky(crossprod(ZtP) + Diagonal(ncol(ZtP), 1), LDL = FALSE)
  	return(list(LV, L))
  }

  cdev <- cmpfun(dev)
  gdev <- cmpfun(grad)
  cL <- cmpfun(cholT)
  list(cdev, gdev, cL)
}
DiACrGG <- cmpfun(MM_BlockDiagGG)

## assume that Lamt is digonal and theta is now log(theta) 
MM_BlockDiagG <- function(y, X, Zt, nZ, s = 1) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)
  Lam <- Diagonal(nrow(Zt), 2)
  mapLam <- rep(1:length(nZ), nZ)
  L0 <- Cholesky(tcrossprod(Zt) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
  w <- rep(1:length(nZ), nZ)
  
  dev <- function(theta) {
    Lam <- Diagonal(nrow(Zt), exp(theta)[mapLam])  ## specific to this model..  
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, -crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	# print(r2)
## Multiply by 2 is part of this.. 
    l1 <- as.numeric(2 * determinant(L)$modulus)
    # print(c(r2, l1, 1/s^2 * r2))
	l1 + (1/s^2 * r2 ) + n * log(2 * pi) + 2 * n*log(s)
  }

  grad <- function(theta) {
    Lam <- Diagonal(nrow(Zt), exp(theta)[mapLam])  ## specific to this model..  
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- backsolve(RX, Xty - crossprod(RXZ, cu), transpose = TRUE)
  	
  	B.hat <- backsolve(RX, cB)
  	u.hat <- solve(L, solve(L, cu - RXZ %*% B.hat, system = 'Lt'), 'P')
  	
  	DL <- solve(L, solve(L, Diagonal(nrow(Zt), 1), system = 'P'), 'L')
    th1 <- rowSums(DL * DL)
    th1 <- nZ - as.numeric(tapply(th1, w, sum))
    th2 <- 1/s^2 * as.numeric(tapply(u.hat * u.hat, w, sum))


	2 * (th1 - th2 )
  }

  
	
  cdev <- cmpfun(dev)
  gdev <- cmpfun(grad)
  list(cdev, gdev)
}


DiACrG <- cmpfun(MM_BlockDiagG)

## nZ just defines teh Zt blocks..
## nZG is a list of pairs.. list(c(1, 4), c(2, 5) etc.. ) that tell me which are correlated
## the pairs are 
MM_Cor <- function(y, X, Zt, nZ, nZG) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)

  ### These are indices for z.. so Z2 is located at sZ[2]:eZ[2]
  ### This will give me an error if I match wrong groups due to diff. lengths
  qZ <- length(nZ)
  nP <- length(nZG)
  eZ <- cumsum(nZ)
  sZ <- c(1, eZ[-qZ] + 1)
  
  ## Do Lam so I can use crossprod Lam, Z
  Lam <- Diagonal(nrow(Zt),  rep(1:qZ, nZ))
  thG <- rep(1, qZ)
  for(i in 1:nP) { ## Lam is lower diag so bigger gets row and smaller gets col
  	j <- nZG[[i]]
  	diag(Lam[sZ[j[2]]:eZ[j[2]], sZ[j[1]]:eZ[j[1]]]) <- i + qZ
    thG[j[2]] <- 1 + i
  }
  x <- Lam@x
  
  L0 <- Cholesky(tcrossprod(crossprod(Lam/(qZ+nP), Zt)) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
  
  ### This decomposition is the Bartlet version in which off diags are not restricted#
  ### and on diags are made positive..  It is simpler w/ derivative but not w/ interpretation and I
  dev <- function(theta) {
  	th <- theta
  	th[1:qZ] <- exp(th[1:qZ])
  	Lam@x <- th[x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	l1 <- as.numeric(2 * determinant(L)$modulus)
    # # print(c(r2, l1, n * (1 + log(2 * pi * r2/n))))

    ## Multiply by 2 is part of this.. 
	c(l1 + n * (1 + log(2 * pi * r2/n)), r2, l1, n * (1 + log(2 * pi * r2/n)))
  }
  
  Cov <- function(theta) {
    d <- diag(exp(theta[1:qZ]))
    for(i in 1:nP) {
      j <- nZG[[i]]
      d[j[2],j[1]] <- theta[qZ + i]
    }
    Sig <- tcrossprod(d)
    s <- sqrt(diag(Sig))
    p <- sapply(nZG, function(j) Sig[j[2],j[1]]/(s[j[1]]*s[j[2]]))
    c(s,p)
  }
  
  cholT <- function(theta) {
  	th <- theta
  	th[1:qZ] <- exp(th[1:qZ])
  	Lam@x <- th[Lam@x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)

  	LV <- Cholesky(crossprod(ZtP) + Diagonal(ncol(ZtP), 1), LDL = FALSE)
  	return(list(LV, L))
  }

  cdev <- cmpfun(dev)
  covP <- cmpfun(Cov)
  cL <- cmpfun(cholT)
  list(cdev, covP, cL)
}
DiaCor <- cmpfun(MM_Cor)



## This one is for cor parameterization.. 
MM_Cor2 <- function(y, X, Zt, nZ, nZG) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)

  ### These are indices for z.. so Z2 is located at sZ[2]:eZ[2]
  ### This will give me an error if I match wrong groups due to diff. lengths
  qZ <- length(nZ)
  nP <- length(nZG)
  eZ <- cumsum(nZ)
  sZ <- c(1, eZ[-qZ] + 1)
  
  ## Do Lam so I can use crossprod Lam, Z
  Lam <- Diagonal(nrow(Zt),  rep(1:qZ, nZ))
  thG <- rep(1, qZ)
  for(i in 1:nP) { ## Lam is lower diag so bigger gets row and smaller gets col
  	j <- nZG[[i]]
  	diag(Lam[sZ[j[2]]:eZ[j[2]], sZ[j[1]]:eZ[j[1]]]) <- i + qZ
    thG[j[2]] <- 1 + i
  }
  w2 <- which(thG > 1)
  x <- Lam@x
  L0 <- Cholesky(tcrossprod(crossprod(Lam/(qZ+nP), Zt)) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
    
    ## This one goes w/ standard parameterization 
  ## with a single correlation, it is not much worse and may work better.. 
  dev <- function(theta) {
  	th <- theta
  	p <- 2*plogis(th[1:nP + qZ])-1
  	th[1:nP + qZ] <- p* exp(th[w2])  	
  	th[1:qZ] <- exp(th[1:qZ]) * (c(1, sqrt(1 - p^2))[thG])

  	
  	Lam@x <- th[x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	l1 <- as.numeric(2 * determinant(L)$modulus)
    # # print(c(r2, l1, n * (1 + log(2 * pi * r2/n))))

    ## Multiply by 2 is part of this.. 
	c(l1 + n * (1 + log(2 * pi * r2/n)), r2, l1, n * (1 + log(2 * pi * r2/n)))
  }
  
  Cov <- function(theta) {
  	p <- 2*plogis(theta[1:nP + qZ])-1
    s <- exp(theta[1:qZ])
    c(s,p)
  }

  cholT <- function(theta) {
  	th <- theta
  	p <- 2*plogis(th[1:nP + qZ])-1

   	th[1:nP + qZ] <- p * exp(th[w2]) 	
  	th[1:qZ] <- exp(th[1:qZ]) * (c(1, sqrt(1 - p^2))[thG])

  	
  	Lam@x <- th[Lam@x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    
  	LV <- Cholesky(crossprod(ZtP) + Diagonal(ncol(ZtP), 1), LDL = FALSE)
  	return(list(LV, L))
  }

  cdev <- cmpfun(dev)
  covP <- cmpfun(Cov)
  cL <- cmpfun(cholT)
  list(cdev, covP, cL)

}
DiaCor2 <- cmpfun(MM_Cor2)







## This one is for cor parameterization.. w/ boundary avoiding..  
MM_Cor3 <- function(y, X, Zt, nZ, nZG) { ##  Z = [Za, Zg, Zb] 
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Zty <- Zt %*% y
  yty <- sum(y * y)	
  n <- length(y)

  ### These are indices for z.. so Z2 is located at sZ[2]:eZ[2]
  ### This will give me an error if I match wrong groups due to diff. lengths
  qZ <- length(nZ)
  nP <- length(nZG)
  eZ <- cumsum(nZ)
  sZ <- c(1, eZ[-qZ] + 1)
  
  ## Do Lam so I can use crossprod Lam, Z
  Lam <- Diagonal(nrow(Zt),  rep(1:qZ, nZ))
  thG <- rep(1, qZ)
  for(i in 1:nP) { ## Lam is lower diag so bigger gets row and smaller gets col
  	j <- nZG[[i]]
  	diag(Lam[sZ[j[2]]:eZ[j[2]], sZ[j[1]]:eZ[j[1]]]) <- i + qZ
    thG[j[2]] <- 1 + i
  }
  w2 <- which(thG > 1)
  x <- Lam@x
  L0 <- Cholesky(tcrossprod(crossprod(Lam/(qZ+nP), Zt)) + Diagonal(nrow(Zt), 1), LDL = FALSE)	
    
    ## This one goes w/ standard parameterization 
  ## with a single correlation, it is not much worse and may work better.. 
  dev <- function(theta) {
  	th <- theta
  	p <- 2*plogis(th[1:nP + qZ])-1
  	th[1:nP + qZ] <- p* exp(th[w2])  	
  	th[1:qZ] <- exp(th[1:qZ]) * (c(1, sqrt(1 - p^2))[thG])
  	
  	Lam@x <- th[x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	l1 <- as.numeric(2 * determinant(L)$modulus)
    # # print(c(r2, l1, n * (1 + log(2 * pi * r2/n))))
     BA <- -sum(dbeta(.5 * (p+1), 2, 2, log=TRUE))
    ## Multiply by 2 is part of this.. 
	c(l1 + n * (1 + log(2 * pi * r2/n)) + 2 * BA, r2, l1, n * (1 + log(2 * pi * r2/n)))
  }

  dev2 <- function(theta) {
  	th <- theta
  	p <- 2*plogis(th[1:nP + qZ])-1
  	th[1:nP + qZ] <- p* exp(th[w2])  	
  	th[1:qZ] <- exp(th[1:qZ]) * (c(1, sqrt(1 - p^2))[thG])
  	
  	Lam@x <- th[x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    RXZ <- solve(L, solve(L, ZtP %*% X, system = 'P'), 'L')
    RX <- chol(XtX - crossprod(RXZ))
  	
  	## SRR specific update
  	cu <- solve(L, solve(L, crossprod(Lam, Zty), system = 'P'), 'L')
    cB <- forwardsolve(t(RX), Xty - crossprod(RXZ, cu))
  	
	r2 <- as.numeric(yty - crossprod(cu) - crossprod(cB))
	l1 <- as.numeric(2 * determinant(L)$modulus)
    # # print(c(r2, l1, n * (1 + log(2 * pi * r2/n))))
     BA <- -sum(dbeta(.5 * (p+1), 2, 2, log=TRUE))
    ## Multiply by 2 is part of this.. 
	c(l1 + n * (1 + log(2 * pi * r2/n)), r2, l1, n * (1 + log(2 * pi * r2/n)))
  }
  
  Cov <- function(theta) {
  	p <- 2*plogis(theta[1:nP + qZ])-1
    s <- exp(theta[1:qZ])
    c(s,p)
  }

  cholT <- function(theta) {
  	th <- theta
  	p <- 2*plogis(th[1:nP + qZ])-1

   	th[1:nP + qZ] <- p * exp(th[w2]) 	
  	th[1:qZ] <- exp(th[1:qZ]) * (c(1, sqrt(1 - p^2))[thG])

  	
  	Lam@x <- th[x]
  	
    ZtP <- crossprod(Lam, Zt)        ### General update
    L <- update(L0, ZtP, mult = 1)
    
  	LV <- Cholesky(crossprod(ZtP) + Diagonal(ncol(ZtP), 1), LDL = FALSE)
  	return(list(LV, L))
  }

  cdev <- cmpfun(dev)
  cdev2 <- cmpfun(dev2)
  covP <- cmpfun(Cov)
  cL <- cmpfun(cholT)
  list(cdev, covP, cdev2, cL)

}
DiaCor3 <- cmpfun(MM_Cor3)

# # ###################              Finally got this workign and it is slower... #######
# system.time(fRL <- lmer(yR ~ D + (1 | Cross) + (1 | Vial), dat, REML=F))
# fRL
# ## i is which theta.. 
# dF <- function(f, theta, i, h = .001) {
  # e1 <- rep(0, length(theta))
  # e1[i] <- h
  
  # t1 <- theta + e1
  # t2 <- theta - e1
  
  
  # 1/(2*h) * (f(t1) - f(t2))
# }

# Z <- cBind(ZCr, ZV)
# nZ <- c(ncol(ZCr), ncol(ZV))

# VC <- as.data.frame(VarCorr(fRL))[c(2,1,3),'sdcor']
# s <- VC[3]
# theta <- log(VC[1:2]/VC[3])
# thF <- exp(theta)^2

# foo <- DiACr(yR, X, t(Z), nZ)
# foo2 <- DiACrGG(yR, X, t(Z), nZ)

# #### Check the values..
# fooI <- nlminb(thF, function(th) foo(th), lower=0)
# thF

# foo(thF)+as.numeric(2*logLik(fRL))
# foo(thF)-foo2[[1]](theta)

# #### Theta value is correct and so everythign else shoud match now also.. 


# ########  Gradients..  
# dL <- list(dF(foo2[[1]], theta+.01, 1), dF(foo2[[1]], theta+.01, 2))
# foo2[[2]](theta+.01) ## they are off in opposite directions.. 
# sapply(dL, function(l) l[1])
# ## fixed effects look fine
# fixef(fRL) - B.foo

# ## All the problems come from the second derivative.. 
# th1 - sapply(dL, function(l) l[3])
# -th2 - sapply(dL, function(l) l[4])  #off by a facorf of 7.65

# plot(unlist(ranef(fRL)[2:1]), u.foo * (exp(theta)[w])) 

# foo <- DiACr(y, X, t(Z), nZ)
# foo2 <- DiACrGG(y, X, t(Z), nZ)

# system.time(f1 <- nlminb(rep(0, length(nZ)), function(th) foo2[[1]](th)))
# system.time(f2 <- nlminb(rep(0, length(nZ)), function(th) foo2[[1]](th), function(th) foo2[[2]](th)))
# system.time(fRL <- lmer(yR ~ D + (1 | Cross) + (1 | Vial), dat, REML=F))


