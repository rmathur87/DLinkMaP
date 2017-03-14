#### This solve generic cross problems where the VC's are all completely diagonal.. 
#### Really just need an X, Zt, and a MapLam
library(compiler)
library(Matrix)
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
## Multiply by 2 is part of this.. 
	as.numeric(2 * determinant(L)$modulus + n * (1 + log(2 * pi * r2/n)))
  }
	
  cdev <- cmpfun(dev)
  cdev
}
DiACr <- cmpfun(MM_BlockDiag)




Rand <- function(Cr1, Cr2) {
  n <- length(Cr1)
  un <- unique(c(Cr1, Cr2))
  p <- length(un)
  
  Z <- matrix(0, n, p, dimnames = list(NULL, Line = un))
  
  for(i in 1:n) Z[i, as.character(c(Cr1[i], Cr2[i]))] <- 1
  
  Z
}

stand <- function(y) (y - mean(y))/sd(y)


rankX <- function(X, Xadd, thres = 1e-3) {
 u <- prcomp(cbind(X[,-1], Xadd))
 X <- cbind(X[,1], u$x[,u$sdev > thres])
 X
}

# # ## Try an example...
# L <- 5000
# Zab <- Diagonal(L)
# Zab[(col(Zab) - row(Zab)) == 1] <- 1
# Zab[L, 1] <- 1
# Zb <- Diagonal(100)[rep(1:100, each = 50),]
# X <- cBind(1, Zab %*% matrix(rnorm(L * 3), L))

# y <- X %*% 1:4 +  Zab %*% rnorm(L, sd = sqrt(2)) + Zb %*% rnorm(100, sd = sqrt(.5)) + 
  # rnorm(L)
# Z <- cBind(Zab, Zb)
# nZ <- c(ncol(Zab), ncol(Zb))


# X[,-1] <- apply(X[,-1], 2, scale)
# f1 <- DiACr(y, X, t(Z), nZ)
# system.time(a2 <- nlminb(c(1, 1), f1, lower = 0))

# f2 <- DiACr(y, X[,1, drop = FALSE], t(Z), nZ)
# system.time(anull <- nlminb(c(1, 1), f2, lower = 0))

# y2 <- scale(y)

# f12 <- DiACr(y2, X, t(Z), nZ)
# system.time(a22 <- nlminb(c(1, 1), f1, lower = 0))

# f22 <- DiACr(y2, X[,1, drop = FALSE], t(Z), nZ)
# system.time(anull2 <- nlminb(c(1, 1), f2, lower = 0))











