args <- commandArgs(trailingOnly = TRUE)

#commDir <- args[1]
commDir <- "Z:/Drosophila/Bayesian/general"
setwd(commDir)
fileArgs <- scan(paste0(commDir, "/parameters.csv"), sep=",", what="")
p <- as.integer(fileArgs[which(fileArgs=='p')+2])
print(paste0("Permutation: ", p))
set.seed(92187)
u <- round(runif(1000) * 2^31)
scriptDir <- fileArgs[4]
print(paste0("Script Dir: ", scriptDir))
outDir <- fileArgs[which(fileArgs=='outDir')+2]
print(paste0("Out Dir: ", outDir))

#.libPaths('/home/sroy/R/x86_64-pc-linux-gnu-library/3.0')
### Do some of the initial function loading just in case.. 
source(paste0(scriptDir, '/FUN.R'))
source(paste0(scriptDir, '/MapFun.R'))
source(paste0(scriptDir, '/gradMM.R'))
### Read in the data...  (Focus on log(y))
source(paste0(scriptDir, '/QTL_Process.R'))
source(paste0(scriptDir, '/MM_Process.R'))
source(paste0(scriptDir, '/QC.R'))
#phenotype <- fileArgs[8]
#print(paste0("Phenotype: ", phenotype))
#fileName <- fileArgs[10]
#print(paste0("File Name: ", fileName))
#dat <- read.table(paste0(commDir, '/', fileName), sep = ',', header = TRUE)
load("Z:/Drosophila/Metabolomic/Davis_Data/finalData.rda")
dat <- final.data

dat$Cross <- factor(dat$cross.type)
#respName <- fileArgs[12]
#print(paste0("Response Name: ", respName))
#y <- dat[,respName]
metaboliteNum <- args[1]
y <- read.table(paste0("Z:/Drosophila/Metabolomic/Davis_Data/artifactResiduals/residual", metaboliteNum, ".csv"), sep=",", header=T)
y <- c(y)$x

ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$female.line), paste0(dat$cross.number, '_', dat$male.line)))


load(paste0(commDir, '/results/metabolites/nullModels/Null', metaboliteNum, '.RData'))

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
LOD <- list()
all.lik <- list()

system.time({
  for(i in 1:nrow(poslist)) {
  #for (i in 5560:5570) {
    #if ((i %% 10) == 0) {
    #  print(i)
    #}
    print(i)
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
    QTL.map <- QTL.F.Map(i, L.V, LineM, LineF, XNull, XLNull, SS.Null, z)
    LOD <- c(LOD, list(QTL.map$out))
    all.lik <- c(all.lik, list(c(i, QTL.map$log.lik)))
    #output:
    #each LOD will be a list of 12 elements
    #LOD[[1]] is the LRT statistic for the Additive Model
    #LOD[[2]] is the LRT statistic for the Dominant Model
    #LOD[[3]] is the LRT statistic for the Full Model
    #LOD[[4]] is the LRT statistic for the Additive-Diet Model
    #LOD[[5]] is the LRT statistic for the Dominant-Diet Model
    #LOD[[6]] is the LRT statistic for the Full-Diet Model
    #LOD[[7]] is the number of cummulative factors for Additive Model
    #LOD[[8]] is the number of cummulative factors for Additive and Dominant Models
    #LOD[[9]] is the number of cummulative factors for Full Model
    #LOD[[10]] is the number of cummulative factors for Additive + Dominant + Full + Add_D
    #LOD[[11]] is the number of cummulative factors for Add + Dominant + Full + Add_D + Dom_D
    #LOD[[12]] is the number of cummulative factors for Add + Dominant + Full + Add_D + Dom_D + Full_D
  }
})

#t(sapply(LOD, LODp)) this finds the p-value associated with each founder strain?
#t(sapply(LOD, LODp)) this finds the p-value associated with each founder strain?
out <- cbind(do.call('rbind', LOD), t(sapply(LOD, LODp)))
colnames(out) <- c("LR - Additive", "LR - Dominant", "LR - Full", "LR - D-Additive", "LR - D-Dominant", "LR - D-Full", 
                   "DF - Additive", "DF - Dominant", "DF - Full", "DF - D-Additive", "DF - D-Dominant", "DF - D-Full",
                   "P-Val - Additivve", "P-Val - Dominant", "P-Val - Full", "P-Val - Main", 
                   "P-Val - Add-D", "P-Val - Dom-D", "P-Val - Full-D", "P-Val - Diet")

like.out <- do.call('rbind', all.lik)
colnames(like.out) <- c("QTL1", "Log Like - Additive", "Log Like - Dominant", "Log Like - Full", 
                        "Log Like - D-Additive", "Log Like - D-Dominant", "Log Like - D-Full")

write.table(out, file = paste0(outDir, '/p-value', metaboliteNum, '.csv'), sep = ',')
write.table(like.out, file=paste0(outDir, '/logLike_', metaboliteNum, '.csv'), sep=',')


