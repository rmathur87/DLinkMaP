######  2014/08/27 the final set up for the TG traits... the LM regression + permutaiton will be easier this
###### way.. 

#.libPaths('/home/sroy/R/x86_64-pc-linux-gnu-library/3.0')
args <- commandArgs(trailingOnly = TRUE)
dir <- "Z:/Drosophila/Bayesian"
source(paste0(dir, '/FUN.R'))
source(paste0(dir, '/MapFun.R'))
source(paste0(dir, '/gradMM.R'))
source(paste0(dir, '/MapFun.R'))
### Read in the data...  (Focus on log(y))
source(paste0(dir, '/QTL_Process.R'))
source(paste0(dir, '/MM_Process.R'))
source(paste0(dir, '/QC.R'))
dat <- read.table(paste0(dir, '/TG_data_updated_5_8_13.txt'), sep = '\t', header = TRUE)
dat <- dat[(dat$Mother != 12338) & (dat$Father != 12338),]
dat$Cross <- factor(dat$Cross)
y <- fLam(.25, dat$TG_concentration)[[1]] * 5
Diet <- (dat$Diet == 'F')  - .5



## For now I could just use cor = 0 and under the null.. 
Zb <- Matrix(model.matrix(~0 + as.factor(cross.number), dat))
ZCr <- Matrix(model.matrix(~ 0 + Cross, dat))
Za <- Matrix(Rand(dat$Mother, dat$Father))
ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$Mother), paste0(dat$cross.number, '_', dat$Father)))
Zg <- Za
for(i in 1:nrow(Zg)) Zg[i, as.character(dat$Mother[i])] <- -1

ZbD <- Zb * Diet
ZM <- Matrix(model.matrix(~0 + as.factor(Month.Crossed), dat))
ZCrD <- ZCr * Diet
ZaD <- Za * Diet
ZBlaD <- ZBla * Diet
ZgD <- Zg * Diet
ZMD <- ZM * Diet
ZV <- Matrix(model.matrix(~0 + as.factor(Vial), dat))
ZP <- Matrix(model.matrix(~0 + as.factor(Plate), dat))
X <- cBind(1, Diet)


Z <- cBind(ZM, Zb, ZBla, ZCr, ZMD, ZbD, ZBlaD , ZCrD, ZV, ZP)
nZ <- c(ncol(ZM), ncol(Zb), ncol(ZBla), ncol(ZCr), ncol(ZMD), ncol(ZbD), ncol(ZBlaD), ncol(ZCrD), ncol(ZV), ncol(ZP))
nZG <- sapply(1:4, function(i) c(i, i + 4), simplify = F)
LineM <- dat$Father
LineF <- dat$Mother

XNull <- X
# fNull <- DiACrGG(y, XNull, t(Z), nZ)
fNull <- DiaCor3(y, XNull, t(Z), nZ, nZG) #This function is completely commented out in the gradMM.R file

# system.time(fNullB <- nlminb(rep(0, length(nZ)), fNull[[1]], fNull[[2]]))
system.time(fNullB <- nlminb(c(rep(0, length(nZ)), rep(0, length(nZG))), fNull[[1]]))
s <- as.data.frame(VarCorr(lmer(y ~ (1|Vial) + (1|Plate), dat)))[3,'sdcor']


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
save(list = c('L.V', 'LineM', 'LineF', 'crossN', 'XNull', 'XLNull', 'SS.Null', 'z'), file = paste0(dir, '/Null_2.RData'))

# LLike: -n/2 * (1+log(2*pi*SS.Null/n)) * 2  
# logLik(lm(z ~ 0+XN)) * 2 is the other way.. 
