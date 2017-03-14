library(DSPRqtl)
library(DSPRqtlDataA)
data(positionlist_wgenetic)
x <- poslist[,4]
bk <- which(poslist[-1,1] != poslist[-nrow(poslist), 1])

##### Read in the data.... 
ReadinDat <- function(f) {
  setwd(f)
  ind <- as.numeric(sapply(strsplit(list.files(), '\\.'), function(x) x[2]))
  nam <- sapply(strsplit(list.files(), '\\.'), function(x) paste(x[1], x[2], 
    x[3], sep = '.'))
  f <- list.files()[order(ind)]
  nam <- nam[order(ind)]


  LOD <- list()

  for(i in 1:length(f)) {
    load(f[i]) 
    assign('Ans', eval(as.name(nam[i])) )
    LOD <- c(LOD, list(Ans))
  }
 

  assign('LOD', do.call('rbind', LOD), envir = .GlobalEnv)
  1
}

## Grps is either the single channel plot or the pair that is subtracted
## 
plotLOD <- function(LOD, Grps, dfNull = 12, main = character(0), 
  x, logp = TRUE) {
	
  nplot <- length(Grps)
  par(mfrow = c(nplot, 1), mar = c(7, 7, 4, 2), mgp = c(4.5, 1, 0))
  ylab <- 'LOD'
  if(logp) ylab = '-log10p'
  for(i in 1:nplot) {
  	y <- LOD[,Grps[[i]]]
    maxy <- max(y + 1, 4)
  	dfy <- LOD[,7+Grps[[i]]] - dfNull
  	
  	if(length(Grps[[i]]) > 1) {
  	   y <- y[,2] - y[,1]
  	   dfy <- dfy[,2] - dfy[,1]
  	   maxy <- max(y + 1, 4)
  	}
  	if(logp) {
  	  y <- -pchisq(2 * log(10) * y, df = dfy, lower.tail = FALSE, log.p=TRUE)/log(10)
  	  maxy <- max(y + 1, 4)
  	}
  	plot(x, y, xlab = "Position (cM)", 
        ylab = ylab, cex.lab = 1.5, type = "n", main = main[i],
        axes = FALSE, ylim = c(0, maxy))
    rect(66.3, -10, 120.3, max(y) + 10, col = "grey90", lty = 0)
    rect(174, -10, 221, max(y) + 10, col = "grey90", lty = 0)
    points(x, y, type = "l", lwd = 2)	
    box()
    lab <- c(0, "66  0", 54, "108  0", 47, 103)
    axis(1, at = c(0, 66.3, 120, 174, 221, 277), labels = lab, tick = FALSE,
       cex.axis = 1)
    axis(2, cex.axis = 1.5)
    mtext("X", line = 2.5, side = 1, at = 33.15, cex = 1.5)
    mtext("2L", line = 2.5, side = 1, at = 93.3, cex = 1.5)
    mtext("2R", line = 2.5, side = 1, at = 147, cex = 1.5)
    mtext("3L", line = 2.5, side = 1, at = 197.5, cex = 1.5)
    mtext("3R", line = 2.5, side = 1, at = 249, cex = 1.5)

  }
}

### This just takes in y values

plotLODGrp <- function(x, y, main = character(0), ylab = character(0), ...) {
	
  par(mar = c(7, 7, 4, 2), mgp = c(4.5, 1, 0))
  matplot(x, y, xlab = "Position (cM)", 
     ylab = ylab, cex.lab = 1.5, type = "n", main = main,
     axes = FALSE, ...)
    rect(66.3, -10, 120.3, max(y) + 10, col = "grey90", lty = 0)
    rect(174, -10, 221, max(y) + 10, col = "grey90", lty = 0)
    for(i in 1:ncol(y)) points(x, y[,i], type = "l", lwd = 2, col = i)	
    box()
    lab <- c(0, "66  0", 54, "108  0", 47, 103)
    axis(1, at = c(0, 66.3, 120, 174, 221, 277), labels = lab, tick = FALSE,
       cex.axis = 1)
    axis(2, cex.axis = 1.5)
    mtext("X", line = 2.5, side = 1, at = 33.15, cex = 1.5)
    mtext("2L", line = 2.5, side = 1, at = 93.3, cex = 1.5)
    mtext("2R", line = 2.5, side = 1, at = 147, cex = 1.5)
    mtext("3L", line = 2.5, side = 1, at = 197.5, cex = 1.5)
    mtext("3R", line = 2.5, side = 1, at = 249, cex = 1.5)
}

Peaks <- function(LOD, df, thres, hthres = 50, logp = TRUE) {
  data(positionlist_wgenetic)
  LOD <- as.matrix(LOD)
  y <- LOD
  if(ncol(LOD) > 1) {
    LOD <- y <- LOD[,2] - LOD[,1]	
  }
  if(logp == TRUE) {
  	y <- -pchisq(2 * log(10) * y, df, lower.tail = FALSE, log.p = TRUE)/log(10)
  }
  chr <- unique(poslist[,1])
  l <- list()
    ### Separate each chr.. 
  for(i in 1:length(chr)) {
    ind <- which(poslist[,1] == chr[i])	
  	indM <- which(y[ind] > thres)
  	w <- numeric(0)
  	if(length(indM) > 0) {
  	  w <- sapply(ind[indM], function(k) {
  	  	rng <- k + -hthres:hthres
  	  	rng <- pmin(pmax(min(ind), rng), max(ind))
  	  	all(y[k] >= y[rng])	
  	  	})
  	}
  	
  	l <- c(l, list(ind[indM[w]]))
  }
  
  
  do.call('c', l)
}

### 2LOD is valid because nearby genomic locations should have similar df.. 

PeaksCI <- function(LOD, w, thres = 2) {
  if(length(w) == 0) return(numeric(0))
  data(positionlist_wgenetic)
  K <- length(w)
  chr <- unique(poslist[,1])
  CIL <- matrix(NA, K, 2)
  for(k in 1:K) {
    LODth <- LOD[w[k]] - thres
    ind <- which(poslist[,1] == poslist[w[k], 1])
    grps <- ind[which(LOD[ind] > LODth)]
    bkpt <- which(diff(grps) != 1)
    grpInd <- findInterval(w[k], c(grps[1], grps[bkpt + 1]))
    CIL[k, 1] <- c(grps[1], grps[bkpt + 1])[grpInd]
    CIL[k, 2] <- c(grps[bkpt], grps[length(grps)])[grpInd]
  }
    
  CIL
}


plotLODPeaks <- function(LOD, Grps, dfNull = 12, main = character(0), 
  x, logp = TRUE, thres = 2.1, LODdrp = rep(2, length(Grps))) {
	
  nplot <- length(Grps)
  par(mfrow = c(nplot, 1), mar = c(7, 7, 4, 2), mgp = c(4.5, 1, 0))
  ylab <- 'LOD'
  if(logp) ylab = '-log10p'
  
  pk <- list()
  pkCI <- list()
  for(i in 1:nplot) {
  	y <- LOD[,Grps[[i]]]
    maxy <- max(y + 1, 4)
  	dfy <- LOD[,7+Grps[[i]]] - dfNull
  	y[is.na(y) | y == 0] <- 0.0001
  	dfy[dfy == 0 | is.na(dfy)] <- 1e10	
  	if(length(Grps[[i]]) > 1) {
  	   y <- y[,2] - y[,1]
  	   dfy <- dfy[,2] - dfy[,1]
  	   maxy <- max(y + 1, 4)
  	}
  	pk <- c(pk, list(Peaks(y, dfy, thres)))
  	pkCI <- c(pkCI, list(PeaksCI(y, pk[[i]], LODdrp[i])))
  	
  	
  	if(logp) {
  	  
  	  y <- -pchisq(2 * log(10) * y, df = dfy, lower.tail = FALSE, log.p=TRUE)/log(10)
  	  maxy <- max(y + 1, 4)
  	}
  	
  	plot(x, y, xlab = "Position (cM)", 
        ylab = ylab, cex.lab = 1.5, type = "n", main = main[i],
        axes = FALSE, ylim = c(0, maxy))
    rect(66.3, -10, 120.3, max(y) + 10, col = "grey90", lty = 0)
    rect(174, -10, 221, max(y) + 10, col = "grey90", lty = 0)
    points(x, y, type = "l", lwd = 2)	
    box()
    lab <- c(0, "66  0", 54, "108  0", 47, 103)
    axis(1, at = c(0, 66.3, 120, 174, 221, 277), labels = lab, tick = FALSE,
       cex.axis = 1)
    axis(2, cex.axis = 1.5)
    if(length(pk[[i]]) > 0) {
      abline(v = x[pkCI[[i]]], lty = 2)
      axis(side = 1, at = x[pk[[i]]], tick = TRUE, labels = FALSE, col = 2)
      pk[[i]] <- cbind(pk[[i]], y[pk[[i]]])
    }
    mtext("X", line = 2.5, side = 1, at = 33.15, cex = 1.5)
    mtext("2L", line = 2.5, side = 1, at = 93.3, cex = 1.5)
    mtext("2R", line = 2.5, side = 1, at = 147, cex = 1.5)
    mtext("3L", line = 2.5, side = 1, at = 197.5, cex = 1.5)
    mtext("3R", line = 2.5, side = 1, at = 249, cex = 1.5)

  }
  
  names(pk) <- names(pkCI) <- main  
  list(pk, pkCI)
}
form <- function(a) {
  format(a, justify = 'left', trim = TRUE, scientific = FALSE)
}

FormatPeaks <- function(PF, PFD) {
  data(positionlist_wgenetic)
  
  nF <- sapply(PF[[1]], nrow)
  nFD <- sapply(PFD[[1]], nrow)
  Eff <- rep(c('F', 'FXD'), c(sum(nF), sum(nFD)))
  Modl <- rep(rep(c('GCA', 'SCA', 'Recip', 'Full'), 2), 
           c(nF, nFD))
  
  Peak <- list(PF, PFD)
  Peak <- do.call('rbind', sapply(Peak, function(x) {
  	cbind(do.call('rbind', x[[1]]), do.call('rbind', x[[2]]))
  }, simplify = FALSE))
           
  PeakL <- paste0(poslist[Peak[,1], 1], ':', 
    form(poslist[Peak[,1], 2]))
  PeakCI <- paste0(poslist[Peak[,1], 1], ':', 
      form(poslist[Peak[,3], 2]), '..', form(poslist[Peak[,4], 2]))
  CI.L <- poslist[Peak[,4], 3] - poslist[Peak[,3], 3]
  
  data.frame(Effect = Eff, Model = Modl, Peak.MB = PeakL, 
   p.value = round(Peak[,2], 3), CI.MB = PeakCI, LenCM = round(CI.L, 3))         
  
}


## Grps are likely founder and founder by diet LOD's.. 
FormatLOD <- function(LOD, grp, dfNull = 12) {
   data(positionlist_wgenetic)
   
   LODy <- sapply(grp, function(x) {
   	 if(length(x) == 1) return(LOD[,x])
   	 LOD[,x[2]] - LOD[,x[1]]
   })
   LODdf <- sapply(grp, function(x) {
   	 if(length(x) == 1) return(LOD[,x+7] - dfNull)
   	 LOD[,x[2]+7] - LOD[,x[1]+7]
   })
   
   LODp <- -pchisq(LODy * 2 * log(10), LODdf, lower.tail = F, log.p = T)/log(10)
   
   out <- data.frame(poslist[,c(1, 2, 4)], LOD = LODy, pVal = LODp, df = LODdf)
   out[1:10,]
   modl <- c('GCA', 'SCA', 'Recip', 'Full')
   colnames(out) <- c('chr', 'Mb', 'Cm', 
     paste0('LOD.', modl), paste0('pVal.', modl), paste0('df.', modl))
   out[,c(1:3, 3 + rep((0:2)*4, 4) + rep(1:4, each = 3))]
}





