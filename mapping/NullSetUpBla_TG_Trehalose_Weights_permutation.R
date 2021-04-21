######  2014/08/27 the final set up for the TG traits... the LM regression + permutaiton will be easier this
###### way..

#.libPaths('/home/sroy/R/x86_64-pc-linux-gnu-library/3.0')
args <- commandArgs(trailingOnly = TRUE)
#dir <- "C:/Users/rmathur/OneDrive - Research Triangle Institute/Documents/GitHub/DLinkMAP/mapping"
dir <- args[1]
setwd(dir)

#fileArgs <- scan(paste0(dir, "/parameters_null.csv"), sep=',', what='')

source("FUN.R")
source("MapFun.R")
source('gradMM.R')
source('QTL_Process.R')
source('MM_Process.R')
source('QC.R')


#Get Values from Parameters File
seed <- args[2] #run seed
print(paste0("Seed: ", seed))
set.seed(seed)
fileName <- args[3] #original phenotype file
print(paste0("FileName: ", fileName))
phenotype <- args[4] #phenotype = 'weight'
print(paste0('Phenotype: ', phenotype))
outDir <- args[5] #output directory
print(paste0('outDir: ', outDir))
numPerm <- args[6] #
print(paste0('Num Permutations: ', numPerm))

if (phenotype == 'metabolite') {
  #metabDir <- fileArgs[which(fileArgs == 'metabDir')+2]
  metabDir <- args[7]
  print(paste0('metabDir: ', metabDir))
} else if (phenotype == 'survival') {
  #survivalType <- fileArgs[which(fileArgs == 'survivalType')+2]
  #print(paste0('survivalType: ', survivalType))
  #survivalTrans <- fileArgs[which(fileArgs == 'survivalTrans')+2]
  #print(paste0('survivalTrans: ', survivalTrans))
  #survivalVarName <- fileArgs[which(fileArgs == 'survivalVarName')+2]
  #print(paste0('survivalVarname: ', surivalVarName))
  print("Survival Analysis is not implemented!!")
}
print("Read Parameters File")

#Read the Dataset
dat <- read.table(fileName, sep = ',', header = TRUE)
print("Read in Dataset")

allPerms = list()
for (p in 1:numPerm) {
    if (phenotype == 'metabolite') {
      metaboliteNum <- args[8]
      #metaboliteNum <- 1
      y <- read.table(paste0(metabDir, "/artifactResiduals/residual", metaboliteNum, ".csv"), sep=",", header=T)
      y <- c(y)$x
    } else {
      dat$y_perm <- sample(dat$y)
      write.table(dat, file = paste0(outDir, '/Data_', phenotype, '_permutation', p, '.csv'), sep = ',', col.names = T, row.names = F, quote = F)
      allPerms = c(allPerms, list(dat$y_perm))
    }
}

#Fixed Effects
dat$Cross <- factor(dat$cross.type)
Diet <- (dat$food.type == '(F) Fat')  - .5
Zb <- Matrix(model.matrix(~0 + as.factor(cross.number), dat))
ZCr <- Matrix(model.matrix(~ 0 + Cross, dat))
Za <- Matrix(Rand(dat$female.line, dat$male.line))
ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$female.line), paste0(dat$cross.number, '_', dat$male.line)))
Zg <- Za
for(i in 1:nrow(Zg)) Zg[i, as.character(dat$female.line[i])] <- -1

ZbD <- Zb * Diet
ZM <- Matrix(model.matrix(~0 + as.factor(month.crossed), dat))
ZCrD <- ZCr * Diet
ZaD <- Za * Diet
ZBlaD <- ZBla * Diet
ZgD <- Zg * Diet
ZMD <- ZM * Diet
ZV <- Matrix(model.matrix(~0 + as.factor(vial), dat))
XNull <- cbind(1, Diet)


if (phenotype == 'tg' | phenotype == 'trehalose') {
  ZP <- Matrix(model.matrix(~0 + as.factor(plate), dat))

  Z <- cBind(ZM, Zb, ZBla, ZCr, ZMD, ZbD, ZBlaD , ZCrD, ZV, ZP)
  nZ <- c(ncol(ZM), ncol(Zb), ncol(ZBla), ncol(ZCr), ncol(ZMD), ncol(ZbD), ncol(ZBlaD), ncol(ZCrD), ncol(ZV), ncol(ZP))

  nZG <- sapply(1:4, function(i) c(i, i + 4), simplify = F)

} else {
  Z <- cBind(ZM, Zb, ZBla, ZCr, ZMD, ZbD, ZBlaD , ZCrD, ZV)
  nZ <- c(ncol(ZM), ncol(Zb), ncol(ZBla), ncol(ZCr), ncol(ZMD), ncol(ZbD), ncol(ZBlaD), ncol(ZCrD), ncol(ZV))
  nZG <- sapply(1:4, function(i) c(i, i + 4), simplify = F)
}
print("Setup Design matrix")

LineM <- dat$male.line
LineF <- dat$female.line

allFnull = list()
for (p in 1:numPerm) {
    if (phenotype == 'Weight') {
      sqtW <- sqrt(as.numeric(dat$number.weighed))

      #fNull <- DiaCor3(sqtW * y, sqtW * XNull, t(sqtW * Z), nZ, nZG)
      thisY <- allPerms[[p]]
      fNull <- DiaCor3(sqtW * thisY, sqtW * XNull, t(sqtW * Z), nZ, nZG)
      allFnull = c(allFnull, fNull)
    } else {
      #fNull <- DiaCor3(y, XNull, t(Z), nZ, nZG)
      thisY <- allPerms[p]
      fNull <- DiaCor3(thisY, XNull, t(Z), nZ, nZG)
    }

    system.time(fNullB <- nlminb(c(rep(0, length(nZ)), rep(0, length(nZG))), fNull[[1]]))


    ##### From the NULL
    L.V <- fNull[[4]](fNullB$par)

    ### Need to residualize w/o a 1.. so avoid using lm or use w/ ~0+
    XLNull <-  solveL(L.V[[1]], XNull)
    XLNull[,2] <- XLNull[,2] - XLNull[,1] * crossprod(XLNull[,1], XLNull[,2])/crossprod(XLNull[,1])
    XN <- XLNull %*% diag(1/sqrt(diag(crossprod(XLNull))))
    z <- as.numeric(solveL(L.V[[1]], thisY))


    SS.0 <- sum(z * z)
    SS.Null <- (1 - sum(colSums(XN * z)^2)/SS.0) * SS.0
    crossN <- dat$cross.number
    # LLike: -n/2 * (1+log(2*pi*SS.Null/n)) * 2
    n <- dim(Z)[1]
    logLike <- -n/2 * (1+log(2*pi*SS.Null/n)) * 2
    print("Computation Complete!! Saving Results!!")

    if (phenotype == 'survival') {
      #save(list= c('L.V', 'LineM', 'LineF', 'crossN', 'XNull', 'XLNull', 'SS.Null', 'z', 'logLike'),
      #     file = paste0(outDir, '/Null_', transform.type, '.RData'))
      print("Survival Analysis is not Implemented!!")
    } else if (phenotype == 'metabolite') {
      save(list= c('L.V', 'LineM', 'LineF', 'crossN', 'XNull', 'XLNull', 'SS.Null', 'z', 'logLike'),
           file = paste0(outDir, '/Null_', metaboliteNum, 'permutation', p, '.RData'))
    } else {
      save(list= c('L.V', 'LineM', 'LineF', 'crossN', 'XNull', 'XLNull', 'SS.Null', 'z', 'logLike'),
           file = paste0(outDir, '/Null_', phenotype, 'permutation', p, '.RData'))
    }
    # LLike: -n/2 * (1+log(2*pi*SS.Null/n)) * 2
    # logLik(lm(z ~ 0+XN)) * 2 is the other way..
}
