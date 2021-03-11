######  2014/08/27 the final set up for the TG traits... the LM regression + permutaiton will be easier this
###### way.. 

#.libPaths('/home/sroy/R/x86_64-pc-linux-gnu-library/3.0')
#args <- commandArgs(trailingOnly = TRUE)
dir <- "~/Mathur_drosophila/general"
source(paste0(dir, '/FUN.R'))
source(paste0(dir, '/MapFun.R'))
source(paste0(dir, '/gradMM.R'))
### Read in the data...  (Focus on log(y))
source(paste0(dir, '/QTL_Process.R'))
source(paste0(dir, '/MM_Process.R'))
source(paste0(dir, '/QC.R'))

#fileName <- "davis_wide_sets_currated5_27_15.txt"
#dat <- read.table(paste0(dir, '/', fileName), sep = ',', header = TRUE)
load("~/Mathur_drosophila/general/finalData.rda")
dat <- final.data

dat$Cross <- factor(dat$cross.type)
#dat$numWeigh <- factor(dat$number.weighed)
#respName <- "trehalose.concentration"
#y <- dat[,respName]
for (metaboliteNum in 219:422) {
  print(paste0("Running Metabolite", metaboliteNum))
  y <- read.table(paste0("~/Mathur_drosophila/general/artifactResiduals/residual", metaboliteNum, ".csv"), sep=",", header=T)
  y <- c(y)$x
  
  Diet <- (dat$food == '(F) Fat')  - .5
  ## For now I could just use cor = 0 and under the null.. 
  Zb <- Matrix(model.matrix(~0 + as.factor(cross.number), dat))
  #Zn <- Matrix(model.matrix(~ 0 + numWeigh, dat))
  ZCr <- Matrix(model.matrix(~ 0 + Cross, dat))
  Za <- Matrix(Rand(dat$female.line, dat$male.line))
  ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$female.line), paste0(dat$cross.number, '_', dat$male.line)))
  Zg <- Za
  for(i in 1:nrow(Zg)) Zg[i, as.character(dat$female.line[i])] <- -1
  
  ZbD <- Zb * Diet
  #ZnD <- Zn * Diet
  ZM <- Matrix(model.matrix(~0 + as.factor(month.crossed), dat))
  ZCrD <- ZCr * Diet
  ZaD <- Za * Diet
  ZBlaD <- ZBla * Diet
  ZgD <- Zg * Diet
  ZMD <- ZM * Diet
  ZV <- Matrix(model.matrix(~0 + as.factor(vial), dat))
  #ZP <- Matrix(model.matrix(~0 + as.factor(plate), dat))
  X <- cBind(1, Diet)
  
  
  #Z <- cBind(ZM, Zb, Zn, ZBla, ZCr, ZMD, ZbD, ZnD, ZBlaD , ZCrD, ZV) #The one for Weights
  Z <- cBind(ZM, Zb, ZBla, ZCr, ZMD, ZbD, ZBlaD , ZCrD, ZV)
  #nZ <- c(ncol(ZM), ncol(Zb), ncol(Zn), ncol(ZBla), ncol(ZCr), ncol(ZMD), ncol(ZbD), ncol(ZnD), ncol(ZBlaD), ncol(ZCrD), ncol(ZV)) #The one for weights
  nZ <- c(ncol(ZM), ncol(Zb), ncol(ZBla), ncol(ZCr), ncol(ZMD), ncol(ZbD), ncol(ZBlaD), ncol(ZCrD), ncol(ZV))
  nZG <- sapply(1:4, function(i) c(i, i + 4), simplify = F)
  LineM <- dat$male.line
  LineF <- dat$female.line
  
  XNull <- X
  # fNull <- DiACrGG(y, XNull, t(Z), nZ)
  fNull <- DiaCor3(y, XNull, t(Z), nZ, nZG)
  
  # system.time(fNullB <- nlminb(rep(0, length(nZ)), fNull[[1]], fNull[[2]]))
  system.time(fNullB <- nlminb(c(rep(0, length(nZ)), rep(0, length(nZG))), fNull[[1]]))
  #s <- as.data.frame(VarCorr(lmer(y ~ (1|vial), dat)))[3,'sdcor']
  
  
  ##### From the NULL 
  L.V <- fNull[[4]](fNullB$par)
  
  ### Need to residualize w/o a 1.. so avoid using lm or use w/ ~0+
  XLNull <-  solveL(L.V[[1]], XNull)
  XLNull[,2] <- XLNull[,2] - XLNull[,1] * crossprod(XLNull[,1], XLNull[,2])/crossprod(XLNull[,1])
  XN <- XLNull %*% diag(1/sqrt(diag(crossprod(XLNull))))
  z <- as.numeric(solveL(L.V[[1]], y))
  
  # This contains extra cor regularization.. 
  # logLik(lm(z ~ 0 + XN)) * 2 - 2 * determinant(L.V[[1]])$modulus + 2 * sum(dbeta(.5 * (fNull[[2]](fNullB$par)[11:14] + 1), 2, 2, log=TRUE))
  
  
  SS.0 <- sum(z * z)
  SS.Null <- (1 - sum(colSums(XN * z)^2)/SS.0) * SS.0
  crossN <- dat$cross.number
  save(list = c('L.V', 'LineM', 'LineF', 'crossN', 'XNull', 'XLNull', 'SS.Null', 'z'), file = paste0(dir, '/results/metabolites/nullModels/Null', metaboliteNum, '.RData'))
  
  # LLike: -n/2 * (1+log(2*pi*SS.Null/n)) * 2  
  # logLik(lm(z ~ 0+XN)) * 2 is the other way.. 
}