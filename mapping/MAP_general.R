#### This file conducts a Mixed Linear Model for a QTL sweep of the Drosophila genome.
#### Inputs: commDir - the main directory where files are stored; 
####         parameters.csv - main parameters for the code


## Read the Arguments from Command Line
args <- commandArgs(trailingOnly = TRUE)
parametersFile = args[1]
parameters = read.table(parametersFile, sep = ',', header = T)

# Script and Seed Parameters
commDir <- parameters[which(parameters[,1]=='commDir'), "Value"] #The directory where dependent scripts are stored
print(paste0("Comm Dir: ", commDir))
setwd(as.character(commDir))
set.seed(parameters[which(parameters[,1]=='seed'), "Value"]) #To ensure that results are replicable
u <- round(runif(1000) * 2^31)

# Phenotype Parameters
outDir <- parameters[which(parameters[,1]=='outDir'), "Value"] #Directory where results will be written to
print(paste0("Out Dir: ", outDir))
phenotype <- parameters[which(parameters[,1]=='phenotype'), "Value"] #Which phenotype is being analyzed - 'metabolite', 'weight', 'trehalose', and 'TG' have been tested
print(paste0("Phenotype: ", phenotype))
weight.sex <- parameters[which(parameters[,1]=='weightSex'), "Value"] #if the 'weight' phenotype is being analyzed, which sex is analyzed - for all other phenotypes, ignore this!
print(paste0('Weight Sex: ', weight.sex))
weight.type <- parameters[which(parameters[,1]=='weightType'), "Value"] # if the 'weight' phenptype is being anayzed, whether the original data or the average by vial is being analyzed - for all other phenotypes, ignore this!!
print(paste0('Weight Type: ', weight.type))
fileName <- parameters[which(parameters[,1]=='fileName'), "Value"] # the filepath and name of the dataset to be analyzed
print(paste0("File Name: ", fileName))
nullDir <- parameters[which(parameters[,1]=='nullDir'), "Value"]
print(paste0('nullDir: ', nullDir))

# Model Parameters
p <- as.integer(parameters[which(parameters[,1]=='p'), "Value"]) #p=0 is no permutation, while p=1 is permutation
print(paste0("Permutation: ", p))
numPermutations <- as.integer(parameters[which(parameters[,1]=='numPermutations'), "Value"]) # Number of simulation to perform for permutation testing
print(paste0('Num Permutations: ', numPermutations))
#Epistatic and Survival Models
epistaticModel <- parameters[which(parameters[,1]=='epistaticModel'), "Value"] # T/F if an epistatic model is run. F is normal mapping with additive, dominant and full effects; T is additive with epistatic effects for all pairs of QTLs
print(paste0("Epistatic Model: ", epistaticModel))
epistaticQTL <- as.integer(parameters[which(parameters[,1]=='epistaticQTL'), "Value"]) # Integer (out of the total number of QTLs tested) for the epistatic model
print(paste0('Epistatic QTL: ', epistaticQTL))
survivalType <- parameters[which(parameters[,1]=='survivalType'), "Value"]
print(paste0('SurivalType: ', survivalType))
survivalTrans <- parameters[which(parameters[,1]=='survivalTrans'), "Value"]
print(paste0('SurvivalTrans: ', survivalTrans))
survivalVarName <- parameters[which(parameters[,1]=='survivalVarName'), "Value"]
print(paste0('survivalVarName: ', survivalVarName))

### Load the dependent R scripts (make sure these are located in the commDir directory)
print(getwd())
source('FUN.R')
source('MapFun_general.R')
source('gradMM.R')
source('QTL_Process.R')
source('MM_Process.R')
source('QC.R')

#Loading the dataset
print("Reading Datafile!!")
dat <- read.table(as.character(fileName), sep = ',', header = TRUE)


if (phenotype == 'survival') {
  print("Survival Data Analysis is not Implemented!!")
  #surv.type <- fileArgs[which(fileArgs == 'survivalType')+2]
  #transform.type <- fileArgs[which(fileArgs == 'survivalTrans')+2]
  #if (transform.type != 'orig') {
  #  load("Transformed.Rda")
  #  surv.var <- fileArgs[which(fileArgs == 'survivalVarName')+2]
  #  y <- get(surv.var)
  
  #  load(paste0(nullDir, "Null_", transform.type, "Trans", ".RData"))
  #} else {
  #  y <- dat$y
  
  #  load(paste0(nullDir, "Null_", transform.type, "Trans", ".RData"))
  #}
  
} else if (phenotype == 'metabolite') {
  print("Running Metabolomic Data Analysis!!")
  metaboliteNum <- as.integer(args[9])
  metabDir <- args[10]
  y <- read.table(paste0(metabDir, "/artifactResiduals/residual", metaboliteNum, ".csv"), sep=",", header=T)
  y <- c(y)$x
  
  load(paste0(nullDir, '/Null_', metaboliteNum, '.RData'))
  
} else if (phenotype == 'weight') {
  print("Loading Null Model Info!!")
  sqrtW <- sqrt(as.numeric(dat$number.weighed))
  y <- dat$y
  
  load(paste0(nullDir, '/Null_', weight.sex, '_', 'Weight_', weight.type, '.RData'))
} else if (phenotype == 'TG' | phenotype == 'trehalose') {
  print("Loading Null Model Info!!")
  y <- dat$y
  
  load(paste0(as.character(nullDir), '/Null2.RData'))
}
print("Loaded in Null Model!!")


dat$Cross <- factor(dat$cross.number) #cross number has to be a factor in the model
#Matrix for the block - mother and father line
ZBla <- Matrix(Rand(paste0(dat$cross.number, '_', dat$female.line), paste0(dat$cross.number, '_', dat$male.line)))
LineM <- dat$male.line
LineF <- dat$female.line


#Permutation Testing - p=0 is no permutation
if(p != 0) {
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
}


######################################           QTL              ######################################
LOD <- list() #List to store the genome LOD (LRT values) and the model DF
all.lik <- list() #List to store the genome log likelihoods - can be used for CI calculations (as done in epistatic model)

print("Starting Genome Scan!!")
system.time({
  for(i in 1:nrow(poslist)) {
    #for(i in 1:10) {
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
    if (as.logical(epistaticModel)) {
      if (i != epistaticQTL) {
        QTL.map <- QTL.F.Map.Inter(i, L.V, LineM, LineF, XNull, XLNull, SS.Null, z, q, sqrtW=sqrtW)
        LOD <- c(LOD, list(QTL.map$out))
        all.lik <- c(all.lik, list(c(i, q, QTL.map$log.lik)))
      } else {
        LOD <- c(LOD, list(rep(NA, 12)))
        all.lik <- c(all.lik, list(c(i, q, rep(NA, 6))))
      }
    } else {
      if(p != 0) {
          LOD_res <- numeric(numPermutations) ## set aside space for results
          lik_res <- numeric(numPermutations)
      }
	    for (perm in 1:numPermutations) {
          QTL.map <- QTL.F.Map(i, L.V, LineM, LineF, XNull, XLNull, SS.Null, z, sqrtW)
          LOD <- c(LOD, list(QTL.map$out))
          all.lik <- c(all.lik, list(c(i, QTL.map$log.lik)))
      }
      QTL.map <- QTL.F.Map(i, L.V, LineM, LineF, XNull, XLNull, SS.Null, z, sqrtW)
      LOD <- c(LOD, list(QTL.map$out))
      all.lik <- c(all.lik, list(c(i, QTL.map$log.lik)))
    }
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
out <- cbind(do.call('rbind', LOD), t(sapply(LOD, LODp))) #Find associated p-values using the LOD scores
like.out <- do.call('rbind', all.lik)

if (epistaticModel) {
  colnames(out) <- c("LR - QTL 1", "LR - QTL 2", "LR - Inter", "LR - D-QTL1", "LR - D-QTL2", "LR - D-Inter", 
                     "DF - QTL 1", "DF - QTL 2", "DF - Inter", "DF - D-QTL1", "DF - D-QTL2", "DF - D-Inter",
                     "P-Val - QTL1", "P-Val - QTL2", "P-Val - Inter", "P-Val - QTL", 
                     "P-Val - QTL1-D", "P-Val - QTL2-D", "P-Val - Inter-D", "P-Val - Diet")
  
  colnames(like.out) <- c('QTL1', 'QTL2', 'Log Like - QTL1', 'Log Like - QTL2', 'Log Like - Inter', 'Log Like - D-QTL1', 'Log Like - D-QTL2', 
                          'Log Like - D-Inter')
  
} else {
  colnames(out) <- c("LR - Additive", "LR - Dominant", "LR - Full", "LR - D-Additive", "LR - D-Dominant", "LR - D-Full", 
                     "DF - Additive", "DF - Dominant", "DF - Full", "DF - D-Additive", "DF - D-Dominant", "DF - D-Full",
                     "P-Val - Additive", "P-Val - Dominant", "P-Val - Full", "P-Val - Main", 
                     "P-Val - Add-D", "P-Val - Dom-D", "P-Val - Full-D", "P-Val - Diet")
  colnames(like.out) <- c("QTL1", "Log Like - Additive", "Log Like - Dominant", "Log Like - Full", 
                          "Log Like - D-Additive", "Log Like - D-Dominant", "Log Like - D-Full")
}

#Write the results
if (epistaticModel) {
  if (phenotype == 'metabolite') {
    write.table(out, file=paste0(outDir, '/p-value_', metaboliteNum, '_', p, '_inter_', q, '.csv'), sep=',')
    write.table(like.out, file=paste0(outDir, '/logLike_', metaboliteNum, '_', p, '_inter_', q, '.csv'), sep=',')
  } else {
    write.table(out, file=paste0(outDir, '/p-value_', p, '_inter_', q, '.csv'), sep = ',')
    write.table(like.out, file=paste0(outDir, '/logLike_', p, '_inter_', q, '.csv'), sep=',')
  }
} else {
  
  if (phenotype == 'metabolite') {
    write.table(out, file = paste0(outDir, '/p-value', metaboliteNum, '.csv'), sep = ',')
    write.table(like.out, file=paste0(outDir, '/logLike_', metaboliteNum, '.csv'), sep=',')
  } else {
    write.table(out, file = paste0(outDir, '/p-value_', weight.sex, '_', weight.type, '.csv'), sep=',')
    write.table(like.out, file=paste0(outDir, '/logLike_', weight.sex, '_', weight.type, '.csv'), sep=',')
  }
}
