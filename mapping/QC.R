## This is to extract the single allelic covariance component
## and the deviance part of it also.. 
library(gplots)
Fnd <- function(pop = c('A', 'B')) {
	if(pop == 'A') {
	  return(matrix(c('A1', 'Canton-S', 0, 0, 
         'A2', 'BOG 1', NA, NA,
         'A3', 'BS 1', 0, 1,
         'A4', 'KSA 2', -1, -1,
         'A5', 'VAG 1', 1, 1,
         'A6', 'Wild 5B', 1, 1,
         'A7', 'tu.7', 1, 1,
         'AB8', 'SAM', NA, NA), ncol = 4, byrow = TRUE) )   
    }
    
    matrix(c('B1', 'BER', NA, NA, 
             'B2', 'CA 1', 1, 1, 
             'B3', 'QI 2', 1, 1, 
             'B4', 'RVC 3', 1, 1, 
             'B5', 'tu.0', 0, 0, 
             'B6', 'tu.1', 0, 1,
             'B7', 'tu.4', 0, 1, 
             'AB8', 'SAM', NA, NA), ncol = 4, byrow = TRUE)
    
}

### This function lets me see how the crosses are spread in terms of 
### each type off possible covariance structure. 

datTab <- function(dat) {
  dat <- sapply(dat, as.character)

  A <- diag(3, nrow(dat), nrow(dat))

  w <- c(1.5, 1.5, 2, 2)
  for(i in 1:(nrow(dat) - 1)) {
    mat <- t(dat[-(1:i), c('F', 'M', 'M', 'F'), drop = FALSE]) == dat[i, c('F', 'M')]
    A[i, -(1:i)] <- colSums(mat * w)
  }

  A[A == 3] <- 'Rep'
  A[A == 4] <- 'Recip'
  A[A == 2] <- 'Line'
  A[A == 1.5] <- 'LineGen'

  outN <- c('Rep', 'Recip', 'Line', 'LineGen')
  out <- table(A)[outN]
  names(out) <- outN
  out[is.na(out)] <- 0
  out
}



coefLME <- function(lmeMod, gInd) {
  Bet <- numeric(0)
  SE <- numeric(0)
  if(gInd[1] != -1) { 
  U1M <- length(fixed.effects(lmeMod)) - length(gInd)
  Bet <- SE <- rep(0, 8)
  Bet[-gInd] <- SE[-gInd] <- NA
  
  ## Assume the first column coefficient is set to 0
  est <- cbind(1, matrix(1/U1M, length(gInd), U1M), rbind(0, diag(length(gInd) - 1)))
  se <- sqrt(diag(est %*% lmeMod$varFix %*% t(est)))
  est <- colSums(t(est) * lmeMod$coefficients$fixed)

  Bet[gInd] <- est
  SE[gInd] <- se
  }
  sig <- lmeMod$sigma
  VC <- as.numeric(sapply(lmeMod$modelStruct$reStruct, function(x) {
  	attr(summary(x), "stdDev")[1]
  }, USE.NAMES = FALSE))
  
  c(Bet, SE, VC * sig, sig)
}

## Grab the Good Indices from X
getGI <- function(X) which((colSums(round(X)) != 0) & (apply(X, 2, max) > .3))

upM <- function(gInd, fixed) {
  Xf <- paste0('X[,', gInd[-1], ']', collapse = " + ")
  upMod <- paste0(fixed, ' + ', Xf)
  upMod
}


### This creates the matrix that will be used to create the 
### Design for dominant model.. 
DomMat <- function(Fn) {
  x1 <- rep(1:Fn, each = Fn)
  x2 <- rep(1:Fn, Fn)
  Dom <- matrix(0, Fn^2, Fn*(Fn + 1)/2)
  col <- 1
  for(i in 1:Fn) {
    for(j in i:Fn) {
      ind <- which(((x1 + x2) == (i + j))  &  
  	    (abs(x1 - x2) == abs(j - i)))
  	  Dom[ind, col] <- 1
  	  col <- col + 1
    }
  }
  Dom
}

## Get pcr vecs using Eigenvalue decomposition incase Lapack doesn't work
pcrEigen <- function(X) {
  X <- sweep(X, 2, colSums(X)/nrow(X))
  V <- eigen(crossprod(X))$vectors
  x <- X %*% V
  
  list(x = x, sdev = apply(x, 2, sd))
}

## Want pcrs cut off is arbitrary but it should be working well.. 
pcrX <- function(X, w = NA, cuto = 1.5e-2) {
  X <- tryCatch(prcomp(X), error = function(e) pcrEigen(X))
  if(is.na(w)) w <- min(which(X$sdev < cuto)) - 1
    
  X$x[,1:w]
}





## This makes a dot plot.. it is probably easier to do it outside function
## but this helps me stay organized.. 
coefPlot <- function(Dat, ind, VCC = NULL, ylab = 'Weight', main, Fndd, rng =1000) {
  Foundr <- Fnd('A')[,2]
  
  dat <- Dat[ind + -rng:rng, 1:4]
  
  
  
  
  x <- as.numeric(Dat[ind, 5:12])
  se <- as.numeric(Dat[ind, 13:20])
  
  
  
  
  
  par(mfrow = c(1, 2))
  plotCI(x, uiw = se, xlab = 'Founder', ylab = 'Effect',
    pch = 18)
  axis(side = 3, at = 1:8, tick = TRUE, labels = Foundr)
  
  x <- poslist[ind + -rng:rng, 'Gpos']
  y <- dat[,4]

  xs <- which.min(x)
  if(xs != 1) {
  	x <- x[-(1:xs)]
  	y <- y[-(1:xs)]
  	
  }
  
  plot(x, y, main = main, xlab = paste0('Position: ', 'cM'),
    ylim = c(0, max(dat[,4]) + 1), type = 'l', ylab = 'LOD')
  leg <- paste0('F', 1:8, ': n = ', as.numeric(Fndd))
  if(!is.null(VCC)) leg <- c(leg, VCC)
  legend('topright', legend = leg)
  
  
}


## Calculate the number of each founder represented in the data
## Entropy is just on lines 
## The number of founders is both on lines 
## The number of crosses is just a binomial mix

Qual <- function(Lines, HMMf, Cut = .95) {
   Probs <- HMMf[, c("A1A1", "A1A2", 
                "A1A3", "A1A4", "A1A5", "A1A6", "A1A7", "A1A8", 
                "A2A2", "A2A3", "A2A4", "A2A5", "A2A6", "A2A7", 
                "A2A8", "A3A3", "A3A4", "A3A5", "A3A6", "A3A7", 
                "A3A8", "A4A4", "A4A5", "A4A6", "A4A7", "A4A8", 
                "A5A5", "A5A6", "A5A7", "A5A8", "A6A6", "A6A7", 
                "A6A8", "A7A7", "A7A8", "A8A8")]
   Hom <- which(substr(colnames(Probs), 1, 2) == substr(colnames(Probs), 3, 4))
   matchLH <- match(Lines, HMMf[,1])             
   ProbsL <- Probs[matchLH, ]
   
   ### Founder Representation
   foundS <- colSums(ProbsL > Cut)
   foundNS <- c(foundS[Hom], 'Het' = sum(foundS[-Hom]), 'Un' = length(matchLH) -sum(foundS))
   
   ### Entropy
    ent <- EntrRIL(ProbsL)
    c(foundNS, ent)
   
}

FSamp <- function(FLines, MLines, HMMf) {
  ProbLine <- HMMf[, 2:9]
  uniFM <- unique(cbind(FLines, MLines))
	
  matchLF <- match(uniFM[,1], HMMf[,1])
  matchLM <- match(uniFM[,2], HMMf[,1])

    Fem <- as.numeric(colSums(t(ProbLine[matchLF, ] > .95) * 1:8))
    Mal <- as.numeric(colSums(t(ProbLine[matchLM, ] > .95) * 1:8))
    drp <- which((Fem == 0) | (Mal == 0) )
    Fem <- Fem[-drp]
    Mal <- Mal[-drp]
    
    O <- matrix(0, 8, 8)
    E <- rep(0, 8)
    for(i in 1:length(Fem))  {
      O[Fem[i], Mal[i]] <- O[Fem[i], Mal[i]] + 1
      E[Mal[i]] <- E[Mal[i]] + 1
      E[Fem[i]] <- E[Fem[i]] + 1
    }
     missE <- which(E == 0) 
    E <- tcrossprod(E[-missE])/sum(E)/2
 
    O <- O[-missE, -missE]
    O <- O + t(O) - diag(diag(O))
    E <- E + t(E) - diag(diag(E))
    O <- O[row(O) >= col(O)]
    E <- E[row(E) >= col(E)]
    TS <- sum(((O - E)^2/E)[E != 0])
    
    TS
}


CrossFound <- function(dat, X) {
  FC <- sapply(unique(sort(dat$CrossNum)), function(i) {
    sum(colSums(round(X[dat$CrossNum == i,])) > 1.9)
  })
  c(sum(colSums(round(X)) > 1.9), FC)
}

Rand <- function(Cr1, Cr2) {
  n <- length(Cr1)
  un <- unique(c(Cr1, Cr2))
  p <- length(un)
  
  Z <- matrix(0, n, p, dimnames = list(NULL, Line = un))
  
  for(i in 1:n) Z[i, as.character(c(Cr1[i], Cr2[i]))] <- 1
  
  Z
}


library(DSPRqtl)
library(DSPRqtlDataA)

Found <- function(dat, i, mcol = 'Mother', pcol = 'Father') {
   data(positionlist_wgenetic)
   objname <- paste("A_", poslist[i, 1], "_", 
                     format(poslist[i, 2], sci = FALSE), sep = "")
   data(list = objname)

    ### Name the data genotypes
    genotypes <- get(objname)
    Z <- Rand(dat[,mcol], dat[,pcol])                
                     
    gind <- match(colnames(Z), genotypes[,1])
    Found <-as.matrix(genotypes[gind, 2:9])
    assign('X', Z %*% Found, envir = .GlobalEnv) 
    Ns <- colSums(round(X))
    ans <- CrossFound(dat, X)            
    rm(list = c(objname), pos = .GlobalEnv) 
    list(ans, Ns = Ns)           
}



EntrRIL <- function(ProbsL) {
  entL <- -(rowSums(ProbsL * log(ProbsL), na.rm = TRUE)/log(8))
  c('Ent' = mean(entL), 'EntSD' = sd(entL))

}


ZCross <- function(Z, Perm) {
  cn <- colnames(Z)[Perm]
  matrix(cn[t(Z) * 1:ncol(Z)], ncol =2, byrow = TRUE)
}


