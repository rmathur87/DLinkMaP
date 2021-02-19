getI <- function(chr, Loc) which(poslist[,1] == chr & poslist[,2] == Loc)

getX <- function(i, LineM, LineF, out = 'Full') {
    objname <- paste("A_", poslist[i, 1], "_", 
                 format(poslist[i, 2], sci = FALSE), sep = "")
    data(list = objname)
    Dom <- DomMat(8)
    ### Name the data genotypes
    genotypes <- get(objname)
    rm(list = objname, pos = .GlobalEnv)
    FoundM <- genotypes[match(LineM, genotypes[,1]), 2:9]
    FoundF <- genotypes[match(LineF, genotypes[,1]), 2:9]
    
    
    
    
	XFull <- as.matrix(FoundM[,rep(1:8, each = 8)] * FoundF[,rep(1:8, 8)]) 

	if(out == 'Full') {
		return(list(Full = XFull, XaddM = as.matrix(FoundM), XaddF = as.matrix(FoundF)))
    }
    if(out == 'SCA') {
        return(list(Full = XFull %*% Dom, XaddM = as.matrix(FoundM), XaddF = as.matrix(FoundF)))
    }
    if(out == 'GCA') {
        return(list(Full = as.matrix(FoundM) + as.matrix(FoundF), 
          XaddM = as.matrix(FoundM), XaddF = as.matrix(FoundF)))
    }
}

datOut <- function(fitM, X, Cross) {
mu <- fixef(fitM)[1] + X %*% fixef(fitM)[-1]

## Out is the data frame of values... 
Out <- data.frame(ranef(fitM)[[2]])
names(Out) <- 'Cross'
Out$CrF <- Out$CrN <- Out$muF <- Out$muN <- Out$Cross

Out$CrF <- ranef(fitM)[[1]][match(paste0('F:', rownames(Out)), rownames(ranef(fitM)[[1]])),]
Out$CrN <- ranef(fitM)[[1]][match(paste0('N:', rownames(Out)), rownames(ranef(fitM)[[1]])),]

dmuF <- unique(data.frame(Cross, mu)[dat$Diet == 'F',])
dmuN <- unique(data.frame(Cross, mu)[dat$Diet == 'N',])

Out$muF <- dmuF[match(rownames(Out), dmuF[,1]),2]
Out$muN <- dmuN[match(rownames(Out), dmuN[,1]),2]

Out2 <- Out[,4:5]
Out2$CrN <- Out$Cross + Out$muN + Out$CrN
Out2$CrF <- Out$Cross + Out$muF + Out$CrF

Out <- data.frame(Avg = .5 * rowSums(Out2), DietEff = Out2$CrF - Out2$CrN,  
                  AvgF = Out2$CrF, AvgN = Out2$CrN)
Out
}


## I have to account for Vial even it that means I can't account for 
## cross by diet as easily..
## so only real VC's are Month and block so remove them then average over..
datOutVial <- function(fitM, dat, REform = ~(1|Vial) + (1|Cross)) {
  yhat <- predict(fitM, dat, REform = REform)

  CrY <- unique(data.frame(dat$Diet:dat$Cross, yhat))
  w <- which(table(as.character(CrY[,1])) > 1)
  for(i in 1:length(w)) {
    ind <- which(CrY[,1] == names(w)[i])
    CrY[ind,2] <- mean(CrY[ind,2])
  }
  dim(unique(CrY))

  Out <- data.frame(ranef(fitM)[['Cross']])
  names(Out) <- 'Cross'
  Out$CrF <- Out$CrN <- Out$Cross

  Out$CrF <- CrY[match(paste0('F:', rownames(Out)), CrY[,1]),2]
  Out$CrN <- CrY[match(paste0('N:', rownames(Out)), CrY[,1]),2]
  
  Out$Avg <- (Out$CrF + Out$CrN)/2
  Out$D.Eff <- Out$CrF - Out$CrN
  
  Out <- Out[,c('Avg', 'D.Eff', 'CrF', 'CrN')]
  colnames(Out)[3:4] <- c('Fat', 'Natural')
  list(Out, yhat)
}


normY <- function(x) {
  qnorm((rank(x) - .5)/length(x))
}

QTLLoc <- function(j, pk) {
  a <- unlist(strsplit(as.character(pk$Peak.MB[j]), ':'))
  i <- which(poslist[,1] == a[1] & poslist[,2] == as.numeric(a[2]))
  eff <- as.character(pk$Effect[j])
  mod <- as.character(pk$Model[j])
  list(i=i, eff=eff, mod=mod, MB=a)
}

XPCR <- function(X, thres =1e-4) {
  Xpcr <- prcomp(X)
  Xout <- Xpcr$x[,Xpcr$sdev > thres]
  list(XMod = Xout, pcr = Xpcr, X=X)
}


## Effect size (rotate back pcr)
B.Eff <- function(XList, FitM) {
  B <- fixef(FitM)
  X <- XList$X
  XMod <- XList$XMod
  
  Xm <- colSums(X)/nrow(X)
  R <- diag(ncol(X) + 1)
  R[1,-1] <- -Xm
  
  RR <- diag(ncol(X)+1)[,0:ncol(XMod)+1]
  RR[-1,-1] <- XList$pcr$rotation[,1:ncol(XMod)]
  BX <- c(R %*% RR %*% B)
  
  list(X = X, B = BX, yhat = cbind(1, X) %*% BX)
}

dat.Eff <- function(dat, XF, Diet, B) {
  CrDat <- unique(data.frame(Cr = dat$Cross, XF))
  X <- CrDat[,-1]
  
  yhatF <- as.matrix(cbind(1, 1, X, X)) %*% as.numeric(B)
  yhatN <- as.matrix(cbind(1, -1, X, -X)) %*% as.numeric(B)
  
  yhat <- matrix(c(yhatF, yhatN), ncol = 2)
  rownames(yhat) <- as.character(CrDat[,1])
  colnames(yhat) <- c('F', 'N')
  yhat
} 

expandPeak <- function(Peak) {
  P <- do.call('rbind', strsplit(as.character(Peak), ':'))
  
  list(chr = P[,1], loc = as.numeric(P[,2]))
}

expandCI <- function(CI) {
  rCI <- do.call('rbind', strsplit(as.character(CI), ':'))
  CI.B <- t(sapply(strsplit(rCI[,2], '\\..'), as.numeric))
  
  list(chr = rCI[,1], loc = CI.B)
}

PinDat <- function(dat1, dat) {
  Pk <- dat1$Peak.MB
  CI <- expandCI(dat$CI.MB)
  P <- expandPeak(Pk)
  a <- rep(NA, length(Pk))
  for(i in 1:length(a)) {
  	a[i] <- any((P[[1]][i] == CI[[1]]) & (P[[2]][i] > CI[[2]][,1]) & 
  	  (P[[2]][i] < CI[[2]][,2]))
  }
  a
}

## Name of file given directory and features..
f.name <- function(d.out, l.qtl){
  paste0(d.out, l.qtl$eff, '_', l.qtl$mod, '_', l.qtl$MB[1], '_', l.qtl$MB[2], '.csv')
}


Anova.MM <- function(fNull, fMain, fInt, X, Xm, Xi) {
  L <- sapply(list(fNull, fMain, fInt), function(x) x[[2]])
  DF <- sapply(list(X, Xm, Xi), function(x) ncol(x))
  p <- pchisq(-diff(L), diff(DF), lower.tail = FALSE, log.p = FALSE)
  out <- cbind(-diff(L), diff(DF),p)
  rownames(out) <- c('Main', 'Int')
  out
}

## returns B For original X, Line, LineD, Cr, CrD, Vial
## s is the fitted vc
## Z is cbind that corresponds to nZ 
## Zout is a list of the componets to extract
## i.e. alphat = tZi Vinv (y - XBhat)
## i.e. 
B.MM <- function(s, Z, nZ, Zout, Zoi = c(3, 7, 4, 8, 9), XList, w = rep(1, nrow(Z))) {
  ### Don't forget the intercept... 
  X <- XList$X
  XMod <- w * cBind(1, XList$XMod)
  
  ## Set up V
  V <- tcrossprod((w * Z) %*% Diagonal(ncol(Z), rep(sqrt(s), nZ))) + Diagonal(nrow(Z), 1)
  R <- Cholesky(V)
  
  ### Get Bhat for the X inputted.. 
  z <- solve(R, y)
  U <- solve(R, XMod)
  
  XVX <- crossprod(XMod, U)
  rx <- chol(.5 * (XVX + t(XVX)))
  XVy <- crossprod(XMod, z)
  ## this is the PCR version.. 
  uhat <- backsolve(rx, backsolve(rx, XVy, transpose = TRUE))
  yhat <- XMod %*% uhat
  
  
  #############  Now output the rest of the stuff...
  ## Start with the random effects.. 
  r <- solve(R, y - yhat)
  alp <- sapply(Zout, function(x) crossprod(w * x, r))
  ## Don't need to change for vial because I only go out to 4
  for(i in 1:length(alp)) alp[[i]] <- alp[[i]] * sqrt(s[Zoi[i]])
  
  #############  Now take care of the fixed effects..  
  ## w doesn't affect anything here because X is transformed first
  ## then weighted..   
  Xm <- colSums(X)/nrow(X)
  R <- diag(ncol(X) + 1)
  R[1,-1] <- -Xm
  
  RR <- diag(ncol(X)+1)[,1:ncol(XMod)]
  RR[-1,-1] <- XList$pcr$rotation[,2:ncol(XMod)-1]
  BX <- c(R %*% RR %*% uhat)

  
  
  list(X = cbind(1, X), B = BX, alp = alp, s = s[Zoi])
}


pred.MM <- function(Bet, Cross, MLine, FLine, Vial, Diet, X, w = rep(1, length(Cross))) {
  ypred.fixed <- cbind(1, Diet, X, X * Diet) %*% Bet$B

  mInd <- match(MLine, rownames(Bet$alp[[1]]))
  fInd <- match(FLine, rownames(Bet$alp[[1]]))  
  CrInd <- match(Cross, substr(rownames(Bet$alp[[3]]), 6, 20))
  VInd <- match(Vial, substr(rownames(Bet$alp[[5]]), 16, 20))
  
  L <- (Bet$alp[[1]][mInd] + Bet$alp[[1]][fInd]) * Bet$s[1] + 
    Diet * (Bet$alp[[2]][mInd] + Bet$alp[[2]][fInd]) * Bet$s[2]
  
  ypred.Line <- ypred.fixed + L
  
  Cr <- Bet$alp[[3]][CrInd] * Bet$s[3] + Diet * Bet$alp[[4]][CrInd] *Bet$s[4] + Bet$alp[[5]][VInd] * Bet$s[5]
  
  ypred.Cross <- ypred.Line + Cr
  
  w * cbind(ypred.fixed, ypred.Line, ypred.Cross)	
}

### Returns the fixed effects differences
### No influence of weight here because I want new values.. 
xbhat.MM <- function(Cross, X, Diet, Bet) {
  Cr <- unique(data.frame(Cross, X))
  X <- as.matrix(Cr[,-1])
  
  yhatF <- as.matrix(cbind(1, 1, X, X) %*% Bet$B)
  yhatN <- as.matrix(cbind(1, -1, X, -X) %*% Bet$B)
  
  yhat <- matrix(c(yhatF, yhatN), ncol = 2)
  rownames(yhat) <- as.character(Cr[,1])
  colnames(yhat) <- c('F', 'N')
  
  data.frame(avg = .5 * rowSums(yhat), D.eff = (yhat[,1] - yhat[,2]), 
    yhat)
}

xbhat.Cross.MM <- function(Cross, MLine, FLine, Vial, Diet, Bet, xb) {
   CrLine <- unique(data.frame(Cross, MLine, FLine))
   mInd <- match(CrLine[[2]], rownames(Bet$alp[[1]]))
   fInd <- match(CrLine[[3]], rownames(Bet$alp[[1]]))  
   CrInd <- match(CrLine[[1]], substr(rownames(Bet$alp[[3]]), 6, 20))

   crhat <- Bet$alp[[3]][CrInd] * Bet$s[3] + Bet$s[4] * cBind(Bet$alp[[4]], -Bet$alp[[4]])[CrInd,]
   Lhat <- (Bet$alp[[1]][mInd] + Bet$alp[[1]][fInd]) * Bet$s[1] + 
              Bet$s[2] * cBind(Bet$alp[[2]][mInd] + Bet$alp[[2]][fInd], 
                            -((Bet$alp[[2]][mInd] + Bet$alp[[2]][fInd])))          
    
   out <- crhat + Lhat
   out <- data.frame(avg = .5 * rowSums(out), D.Eff = out[,1] - out[,2], F = out[,1], N = out[,2])
   
   out <- out + xb
   rownames(out) <- rownames(xb)
   out
}


reducePk <- function(pk, cut = .001) {
  wMain <- which(pk$pMain < cut)
  wInt <- which(pk$pInt < cut)
  
  mPMain <- min(pk$pMain[pk$pMain > cut])
  mPInt <- min(pk$pInt[pk$pInt > cut])
  
  list(Main = wMain, Int = wInt, c(mPMain, mPInt))
}


makeDietD <- function(wInt, xb) {
  cl <- grep('D.eff', colnames(xb[[1]]))
  l <- list()
  for(i in 1:length(cl)) {
  	out <- sapply(xb[wInt], function(x) x[,cl[i]])
  	colnames(out) <- names(xb[wInt])
  	rownames(out) <- rownames(xb[[1]])
  	l <- c(l, list(out))
  }
  nam <- c('fixed', 'cross')
  names(l) <- nam[1:length(l)]
  l
}

xb.out <- function(dout, w, xb) {  
  ## each full and int gets a value..
  ## label Full.Name.csv or Int.Name.csv
  ## Right now just do the xb but add crosses in a moment
  
  ## Main
  sapply(w[['Main']], function(i) {
    write.csv(xb[[i]], file = paste0(dout, 'Main.', names(xb)[i], '.csv'))
  })

  sapply(w[['Int']], function(i) {
    write.csv(xb[[i]], file = paste0(dout, 'Int.', names(xb)[i], '.csv'))
  })
  
  1
} 



























