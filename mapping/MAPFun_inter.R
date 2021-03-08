fileArgs <- scan(paste0(dir, '/parameters_varExp_inter.csv'), sep=',', what='')
scriptDir <- fileArgs[which(fileArgs == 'scriptDir')+2]
#dir <- "Z:/Drosophila/Bayesian"
source(paste0(scriptDir, '/QTL_Process.R'))
source(paste0(scriptDir, '/MM_Process.R'))
source(paste0(scriptDir, '/QC.R'))



library(DSPRqtl)
library(DSPRqtlDataA)
data(positionlist_wgenetic)
library(lme4)
### This grabs all 3 froms of X from the DSPR hap positions.. 
makeX <- function(i, LineM, LineF) {
  
  objname <- paste("A_", poslist[i, 1], "_", 
                   format(poslist[i, 2], sci = FALSE), sep = "")
  data(list = objname)
  Dom <- DomMat(8)
  ### Name the data genotypes
  genotypes <- get(objname)
  #genotypes <- hmmA[which((hmmA[,1] == poslist[i,1]) & (hmmA[,2] == poslist[i,2])),]
  rm(list = objname, pos = .GlobalEnv)
  FoundM <- genotypes[match(LineM, genotypes[,1]), 2:9] # matches the line of the maternal with this line
  FoundF <- genotypes[match(LineF, genotypes[,1]), 2:9] # matches the line of the paternal with this line  
  
  XFull <- as.matrix(FoundM[,rep(1:8, each = 8)] * FoundF[,rep(1:8, 8)]) 
  XDom <- XFull %*% Dom
  XAdd <- as.matrix(FoundM) + as.matrix(FoundF)
  X <- cBind(XAdd, XDom, XFull)
  nZ <- c(8, ncol(XDom), ncol(XFull))
  list(XA = XAdd, XD = XDom, XF = XFull, FoundM = FoundM, FoundF = FoundF)  
}

### This provides indices matches to turn Full into dom into additive...
### the math formulas are pij_F = pi * pj, pij_D = (pi pj + p1j p2i) if i neq j, and pij_a = pi + pj
## that is sloppy but it should be clear.. 
LF_Mod <- function(L) {
  fi <- rep(1:L, each = L)
  fj <- rep(1:L, L)
  
  
  di <- rep(1:L, L:1)
  dj <- unlist(sapply(1:L, function(i) i:L))
  
  fd <- match(paste0(pmin(fi,fj), '-', pmax(fi,fj)), paste0(di, '-', dj))
  
  list('fD' = fd, 'DAi' = di, 'DAj' = dj)
}

### the box cox transformation for lam along with the modificaiton to LLike.. 
fLam <- function(lam, y) {
  if(lam == 0) return(list('y' = log(y), 'dlgp' = -sum(log(y))))
  return(list('y' = (y^lam - 1)/lam, 'dlgp' = sum((lam - 1) * log(y))) )
}

## Box cox inverse transformation
fILam <- function(lam, y) {
  if(lam == 0) return(exp(y))
  (lam * y + 1)^(1/lam)
}

## projection matrix for dense matrices using QR decompositon.. 
## rank gets determined automatically in the decomposition but I have to be
## careful 
px <- function(x) {  ## QQ^T = PX 
  Q <- qr(x)    
  qr.Q(Q)[,1:Q$rank]
}

## for sprase matrix, teh QR is a different object and so this function is
## required.. the tol should be changed as needed.. 
px.S <- function(x) { ## Sparse matrix is annoying
  Q <- qr(x)	
  R0 <- drop0(qrR(Q), tol = 1e-10) 
  ## columns of R are permuted but rows that are 0 correspond to bad columns of Q
  w <- which(rowSums(R0) != 0)
  qr.Q(Q)[,w,drop=F]
}

## the main gola of the projection matrix: densee
pyhatx <- function(x, y) {
  Q <- px(x)
  DF <- ncol(Q)
  list(as.numeric(Q %*% crossprod(Q, y)), DF)
}

## the sparse version
pyhatSx <- function(x, y) {
  Q <- px.S(x)
  DF <- ncol(Q)
  list(as.numeric(Q %*% crossprod(Q, y)), DF)
}

## S is used here to make sure that crossprod(X) gives I
## X'S are sclaed so that I maintain interpretable coefficients..
## column space is unaffected so the rest is irrelevant..
## this combines each matrix into one large one as I proceed through them
## this would be useful for fitting the full model
fX <- function(XL, X, S = diag(crossprod(X)), thres = 1e-3) {
  if(length(dim(X)) == 2) n <- nrow(X)
  if(length(dim(X)) == 0) n <- length(X)
  U <- X %*% diag(1/sqrt(S))
  S <- S
  Xout <- list()
  
  Min <- apply(U, 2, min)
  Max <- apply(U, 2, max)
  
  for(i in 1:length(XL)) {
    X <- svd(XL[[i]] - U %*% crossprod(U, XL[[i]]))
    w <- which(X$d > thres)
    Min <- c(Min, apply(X$u[,w], 2, min))
    Max <- c(Max, apply(X$u[,w], 2, max))
    U <- cbind(U, X$u[,w])
    
    Xout <- c(Xout, list(as(U %*% Diagonal(ncol(U), pmin(50, 1/(Max-Min))), 'matrix')))
  }
  
  Xout
}


## S is used here to make sure that crossprod(X) gives I
## X'S are sclaed so that I maintain interpretable coefficients..
## column space is unaffected so the rest is irrelevant..
## This funciton just adds in the new X's that are used..
## in orthogonal X this is simpler.. and for LRT's it is all that is needed.. 
fX.R2 <- function(XL, X, S = diag(crossprod(X)), thres = 1e-3) {
  if(length(dim(X)) == 2) n <- nrow(X)
  if(length(dim(X)) == 0) n <- length(X)
  U <- X %*% diag(1/sqrt(S))
  S <- S
  Xout <- list()
  
  Min <- apply(U, 2, min)
  Max <- apply(U, 2, max)
  
  for(i in 1:length(XL)) {
    X <- svd(XL[[i]] - U %*% crossprod(U, XL[[i]]))
    w <- which(X$d > thres)
    Min <- c(Min, wm <- apply(X$u[,w], 2, min))
    Max <- c(Max, wM <- apply(X$u[,w], 2, max))
    U <- cbind(U, X$u[,w])
    #%*% Diagonal(length(w), pmin(50, 1/(wM - wm))
    Xout <- c(Xout, list(as(X$u[,w], 'matrix')))
  }
  
  Xout
}
## changed svd to prcomp because it seems to work better.. 
## threshold may make more sense at this level also.. 
##make sure X is orthogonal
fX.R <- function(XL, X, S = diag(crossprod(X)), thres = 1e-6) {
  if(length(dim(X)) == 2) n <- nrow(X)
  if(length(dim(X)) == 0) n <- length(X)
  U <- X %*% diag(1/sqrt(S))
  S <- S
  Xout <- list()
  for(i in 1:length(XL)) {
    X <- prcomp(XL[[i]] - U %*% crossprod(U, XL[[i]]), center = FALSE)
    w <- which(X$sdev > thres)
    u <- X$x[,w,drop=F] %*% diag(1/sqrt(colSums(X$x[,w,drop=F]^2)))
    U <- cbind(U, u)
    #%*% Diagonal(length(w), pmin(50, 1/(wM - wm))
    Xout <- c(Xout, list(as(u, 'matrix')))
  }
  Xout
}

## the solver for sparse matrices.. MAKE SURE TO NOT USE LDL!!
solveL <- function(L, x){
  as(solve(L, solve(L, x, system = 'P'), 'L'), 'matrix')
}

## teh LOD to p converter for output from QTL.F.Map
LODp <- function(out) {
  dL <- c(out[1], diff(out[1:6]))
  dL <- c(dL[1:3], out[3], dL[4:6], out[6]-out[3])
  
  dfL <- c(out[7], diff(out[7:12]))
  dfL <- c(dfL[1:3], out[9], dfL[4:6], out[12]-out[9])
  -pchisq(2*dL, dfL, lower.tail = FALSE, log.p=T)/log(10)
}

### The greatly simplified function with extra components because they won't be defined in 
## the workspace.. 
QTL.F.Map <- function(QTLI, L.V, LineM, LineF, XNull, XLNull, SS.Null, z, interQTL, sqrtW=1) {
  #QTLI is the position of the qtl
  #L.V is something about the null model
  #LineM is the line of the father
  #LineF is the line of the mother
  #XNull is the design matrix under the null
  #XLNull utilizes the sparse matrix solver
  #SS.Null is the sum of squares under the null
  #z utilizes the sparse matrix solver
  
  #if (QTLI == interQTL) {
    #return(rep(NA, 16))
  #}
  
  ### Set Up X
  n <- length(z)
  XQTL1 <- makeX(QTLI, LineM, LineF) # grabs the additive, dominance and founder  haplotype values
  XQTL2 <- makeX(interQTL, LineM, LineF) #grabs the additive, dominance and founder haplotype values for the interacting QTL
  interact <- model.matrix(~0+XQTL1[[1]] * XQTL2[[1]])[,-(1:16)] #Interaction Effect. Taking out the first 16 since these are the main effects.
  if (QTLI == interQTL) {
    XQTL <- list(XQTL1[[1]])
    Xq <- c(XQTL[1:3], sapply(XQTL[1:3], function(x) list(XNull[,2] * x))) #this QTL x Diet Interactions (not including the interaction term for this)
    names(Xq) <- c("QTL1_A", "QTL1_AD")
    XL <- sapply(Xq, function(x) (solveL(L.V[[1]], x) * sqrtW)) #transforming qtl part of x
    #XL <- sapply(Xq, function(x) solveL(L.V[[1]], x)) #transforming qtl part of x
    XXR <- tryCatch(fX.R(XL, XLNull), error = function(e) {print('LaPack!!');NA}) #PCA part
    if(is.na(XXR[1])) return(list(out=rep(NA, 12), lik=rep(NA, 6)))
  } else {
    XQTL <- list(XQTL1[[1]], XQTL2[[1]], interact)
    Xq <- c(XQTL[1:3], sapply(XQTL[1:3], function(x) list(XNull[,2] * x))) #this QTL x Diet Interactions (not including the interaction term for this)
    names(Xq) <- c("QTL1_A", "QTL2_A", "QTL12_A", "QTL1_AD", "QTL2_AD", "QTL12_AD")
    XL <- sapply(Xq, function(x) (solveL(L.V[[1]], x) * sqrtW)) #transforming qtl part of x
    #XL <- sapply(Xq, function(x) solveL(L.V[[1]], x)) #transforming qtl part of x
    XXR <- tryCatch(fX.R(XL, XLNull), error = function(e) {print('LaPack!!');NA}) #PCA part
    if(is.na(XXR[1])) return(list(out=rep(NA, 12), lik=rep(NA, 6)))
  }
  
  ### Compute LRT.. -n/2(1 + log(2 * pi * SS/n))
  ### since LRT-diff is all that is needed.. I subtract from the null
  ### So LRTa - LRTn = -n * log(SSa/SSn) = -n * log(1-SS(LOFa)/SSn) due to orthgonal nature..
  Ly <- sapply(XXR, function(x) tryCatch(sum(colSums(x * z)^2), error = function(e) {NA}))
  
  ### LRT and the df.. 
  L <- -n/2 * log1p(-cumsum(Ly)/SS.Null)
  #L <- c(L, (L[3]-L[2]), (L[6]-L[5]))
  out <- c(L, cumsum(sapply(XXR, ncol)))
  #out <- c(out, (out[11]-out[10]), (out[14]-out[13]))
  
  SSEm <- SS.Null - cumsum(Ly)
  likelihood <- (-n/2) * (1+log(2*pi*(SSEm/n)))
  
  list(out=out, lik=likelihood) #LRT test statistic
}

### Function used to shuffle LineM and Line from Father
### This can take care of the whole data or loop over blocks.. 
BlS <- function(LM, LF) {    
  U <- unique(c(LM, LF))
  LMi <- match(LM, U)
  LFi <- match(LF, U)
  Up <- sample(U)
  LM <- Up[LMi]
  LF <- Up[LFi]
  list(LM, LF)
}