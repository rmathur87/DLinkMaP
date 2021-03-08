### This file conducts a Epistatic Mixed Linear Model for a pairwise sweep of the Drosophila genome.
### Inputs: commDir - the main directory where files are stored.
###         parameters.csv - main parameters for the code


args <- commandArgs(trailingOnly = TRUE)

commDir <- args[1]
#commDir <- "Z:/Drosophila/Bayesian/general/"
setwd(commDir)

fileArgs <- scan("parameters_inter.csv", sep=",", what="")
p <- as.integer(fileArgs[which(fileArgs == 'p')+2])
print(paste0("Permutation: ", p))
set.seed(92187)
u <- round(runif(1000) * 2^31)

q <- as.integer(args[2])
#q <- 2
print(paste0("Interacting QTL: ", q))

#.libPaths('/home/sroy/R/x86_64-pc-linux-gnu-library/3.0')
### Do some of the initial function loading just in case.. 
source('FUN.R')
source('MAPFun_inter.R')
source('gradMM.R')
### Read in the data...  (Focus on log(y))
source('QTL_Process.R')
source('MM_Process.R')
source('QC.R')
phenotype <- fileArgs[which(fileArgs == 'phenotype')+2]
print(paste0('Phenotype: ', phenotype))
fileName <- fileArgs[which(fileArgs == 'fileName')+2]
print(paste0('File Name: ', fileName))

#Loading the dataset
print("Reading Datafile!!")
dat <- read.table(fileName, sep = ',', header = TRUE)

#Defining the Response and Loading the Null Model (Null Model is created by different script)
nullDir <- fileArgs[which(fileArgs == 'nullDir')+2]
weight.sex <- fileArgs[which(fileArgs == 'weightSex')+2]
weight.type <- fileArgs[which(fileArgs == 'weightType')+2]
if (phenotype == 'metabolite') {
  print('Running Metabolomic Data Analysis!!')
  metaboliteNum <- args[2]
  metabDir <- fileArgs[which(fileArgs == 'metabDir')+2]
  y <- read.table(paste0(metabDir, '/artifactResiduals/residual', metaboliteNum, '.csv'), sep=',', header=T)
  y <- c(y)$x
  
  load(paste0(nullDir, '/Null_', metaboliteNum, '.RData'))
  
} else if (phenotype == 'weight') {
  print("Loading Null Model Info!!")
  sqrtW <- sqrt(as.numeric(dat$number.weighed))
  y <- dat$y
  
  load(paste0(nullDir, '/Null_', weight.sex, '_', 'Weight_', weight.type, '.RData'))
} else {
  print('Loading Null Model Info!!')
  y <- dat$y
  
  load(paste0(nullDir, '/Null2.RData'))
}


dat$Cross <- factor(dat$cross.number) #cross number has to be a factor in the model
#Matrix for the block - mother and father line
ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$female.line), paste0(dat$cross.number, '_', dat$male.line)))


if(p != 0) { #p=0 corresponds to no permutation.. 
  set.seed(u[p])
  LM <- LineM
  LF <- LineF
  cn <- unique(dat$cross.number)
  for(i in 1:length(cn)) {
    w <- which(dat$cross.number == cn[i])
    LM.O <- BlS(LM[w], LF[w])
    LM[w] <- LM.O[[1]]
    LF[w] <- LM.O[[2]]
  }
  LineM <- LM
  LineF <- LF
# ZBla2 <- Matrix(Rand(paste0(dat$cross.number, '_', LM), paste0(dat$cross.number, '_', LF)))
# all(tcrossprod(ZBla2) == tcrossprod(ZBla))
}


######################################           QTL              ######################################
LOD <- list() #List to store the genome LOD (LRT values) and the model DF
all.lik <- list() #List to store the geome log likelihoods - can be used for CI calculations (as done in epistatic model)

print('Starting Genome Scan!!')
system.time({
  #for(i in q:nrow(poslist)) {
  for (i in 1:10) {
    if ((i %% 10) == 0) {
      print(i)
    }
    #Input:
    #i is the position in the genome
    #LineM is the father information. This is defined in the NullSetUpBla.R file
    #LineF is the mother information. This is defined in the NullSetUpBla.R file
    #XNull is the design matrix under the null - only diet
    #XLNull utilizes the sparse matrix solver
    #SS.Null is the sum of squares under the null
    #z is the cholesky's decomposition under the null - only random effects
    #L.V is the cholesky decomposition of V
    #q is the QTL position that will be tested for interaction
    if (i != q) {
      QTL.map <- QTL.F.Map(i, L.V, LineM, LineF, XNull, XLNull, SS.Null, z, q, sqrtW=sqrtW)
      LOD <- c(LOD, list(QTL.map$out))
      all.lik <- c(all.lik, list(c(i, q, QTL.map$log.lik)))
    } else {
      LOD <- c(LOD, list(rep(NA, 12)))
      all.lik <- c(all.lik, list(c(i, q, rep(NA, 6))))
    }
    #output:
    #each LOD will be a list of 12 elements
    #LOD[[1]] is the LRT statistic for QTL1
    #LOD[[2]] is the LRT statistic for QTL2
    #LOD[[3]] is the LRT statistic for QTL12
    #LOD[[4]] is the LRT statistic for QTL1_D
    #LOD[[5]] is the LRT statistic for QTL2_D
    #LOD[[6]] is the LRT statistic for QTL12_D
    #LOD[[7]] is the number of cummulative factors for QTL1
    #LOD[[8]] is the number of cummulative factors for QTL1 + QTL2
    #LOD[[9]] is the number of cummulative factors for QTL1 + QTL2 + QTL12
    #LOD[[10]] is the number of cummulative factors for QTL1 + QTL2 + QTL12 + QTL1_D
    #LOD[[11]] is the number of cummulative factors for QTL1 + QTL2 + QTL12 + QTL1_D + QTL2_D
    #LOD[[12]] is the number of cummulative factors for QTL1 + QTL2 + QTL12 + QTL1_D + QTL2_D + QTL12_D
  }
})

#t(sapply(LOD, LODp)) this finds the p-value associated with each founder strain?
out <- cbind(do.call('rbind', LOD), t(sapply(LOD, LODp)))
colnames(out) <- c("LR - QTL 1", "LR - QTL 2", "LR - Inter", "LR - D-QTL1", "LR - D-QTL2", "LR - D-Inter", 
                   "DF - QTL 1", "DF - QTL 2", "DF - Inter", "DF - D-QTL1", "DF - D-QTL2", "DF - D-Inter",
                   "P-Val - QTL1", "P-Val - QTL2", "P-Val - Inter", "P-Val - QTL", 
                   "P-Val - QTL1-D", "P-Val - QTL2-D", "P-Val - Inter-D", "P-Val - Diet")

like.out <- do.call('rbind', all.lik)
colnames(like.out) <- c('QTL1', 'QTL2', 'Log Like - QTL1', 'Log Like - QTL2', 'Log Like - Inter', 'Log Like - D-QTL1', 'Log Like - D-QTL2', 
                        'Log Like - D-Inter')


#Write the Results
outDir <- fileArgs[which(fileArgs == 'outDir')+2]
if (phenotype == 'metabolite') {
  write.table(out, file=paste0(outDir, '/p-value_', metaboliteNum, '_', p, '_inter_', q, '.csv'), sep=',')
  write.table(like.out, file=paste0(outDir, '/logLike_', metaboliteNum, '_', p, '_inter_', q, '.csv'), sep=',')
} else {
  write.table(out, file=paste0(outDir, '/p-value_', p, '_inter_', q, '.csv'), sep = ',')
  write.table(like.out, file=paste0(outDir, '/logLike_', p, '_inter_', q, '.csv'), sep=',')
}



