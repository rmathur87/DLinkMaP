library(DSPRqtl)
library(ggplot2)
library(hash)
source('~/Drosophila/Bayesian/CI/CI/fun.R')

manhattanPlot <- function(chr.vec, bp.vec, pVal.vec, modelType, dataType, outDir, peakPos=NA, thres=6) {
  P <- pVal.vec[!is.na(pVal.vec)]
  CHR <- chr.vec[!is.na(pVal.vec)]
  BP <- bp.vec[!is.na(pVal.vec)]/1e6
  theData <- data.frame(xAxis=BP, yAxis=P, chr=CHR)
  
  pointsColors <- rep(1, length(P))
  if (!is.na(peakPos)) {
    for (aPeak in 1:nrow(peakPos)) {
      peakInd <- which((CHR == as.character(peakPos[aPeak,1])) & (BP == (peakPos[aPeak,2]/1e6)))
      pointsColors[peakInd] <- 0
    }
  }
  theData$pointColors <- pointsColors
  
  #g1 <- ggplot(theData) + geom_point(aes(x=xAxis,y=yAxis,group=chr)) + 
  #  facet_grid(.~chr,scales="free") + geom_hline(yintercept = 6, col='red') +
  #  labs(title=paste0(dataType, " - ", modelType),x="QTL Position (mB)",y="-log10(p)")
  g1 <- ggplot(theData) + geom_point(aes(x=xAxis,y=yAxis,group=chr,colour=factor(pointColors))) + 
    facet_grid(.~chr,scales="free") + geom_hline(yintercept = thres, col='red') +
    labs(title=paste0(dataType, " - ", modelType),x="QTL Position (Mb)",y="-log10(p)")
  g1 <- g1 + theme(legend.position="none")
  ggsave(filename=paste0(outDir, "/manhattan_", dataType, "_", modelType, "2.pdf"), plot=g1)
}

pk_find <- function(p, thres, model, cmRange = 2) {
  
  pk <- rep(0, nrow(poslist))
  
  for(i in 1:nrow(poslist)) {
    ## +- 2cm
    ind <- which(abs(poslist[,4] - poslist[i, 4]) < cmRange &
                   poslist[,1] == poslist[i, 1])
    if(is.na(p[i])) {
      pk[i] <- 0
    } else {
      pk[i] <- all(p[i] >= p[ind], na.rm = TRUE) & p[i] > thres
    }
    #print(i) 
  }
  w <- which(pk == 1)
  #plot(poslist[,4], p, pch = 20, xlab="Pos (cM)", ylab="-log(p)", main=paste0("Manhattan Plot - ", model, " (thres=", thres, ")")) # ugly LOD plot
  #points(poslist[w, 4], p[w], pch = 20, col = 2)
  
  return(w)
}

theDir <- "~/Drosophila/Bayesian/general/results/metabolites/updated/"
writePeaks.ind <- F
eachManhattan <- F
findCIs <- F

CHR <- poslist[,1] #1="2L"; 2="2R"; 3="3L"; 4="3R"; 5="X"
BP <- as.numeric(poslist[,2])
Gpos <- as.numeric(poslist[,3])

metab.names <- colnames(read.csv("~/Drosophila/Metabolomic/Davis_Data/finalData.csv", sep=',', header=T))[19:440]

#For overlaying manhattan plots - use matrices to store the data
P.Val...Additivve.matrix <- NULL
P.Val...Dominant.matrix <- NULL
P.Val...Full.matrix <- NULL
P.Val...Main.matrix <- NULL
P.Val...Add.D.matrix <- NULL
P.Val...Dom.D.matrix <- NULL
P.Val...Full.D.matrix <- NULL
P.Val...Diet.matrix <- NULL

#For comparing peaks between metabolites - use matrices for main and diet peaks
#allMainPeaks.file <- paste0(theDir, 'plots/allMainPeaks2.csv')
#write('Model Type,Metab Ind,Metab Name,QTL Ind,QTL chr,QTL Pos(Mb),QTL CI lower,QTL CI upper,QTL CI length,p-value', file=allMainPeaks.file, append=FALSE)
#allDietPeaks.file <- paste0(theDir, 'plots/allDietPeaks2.csv')
#write('Model Type,Metab Ind,Metab Name,QTL Ind,QTL chr,QTL Pos(Mb),QTL CI lower,QTL CI upper,QTL CI length,p-value', file=allDietPeaks.file, append=FALSE)

metabPeaks.file <- paste0(theDir, 'plots/metabPeaks.txt')
write('Metabolite\tModel\tNum Peaks\tQTL Indices', file=metabPeaks.file, append=F)

for (i in 1:length(metab.names)) {
#for (i in 1:10) {
  print(paste0("Running Metabolite ", i))
  results <- tryCatch({
    read.csv(paste0(theDir, "genomeScan/p-value", i, ".csv"), header=T, sep=",")},
    error=function(e) {results = NULL})
  
  if (!is.null(results)) {
    diffPVals <- names(results)[grepl("P.Val", names(results))]
    fileNames <- c("Additive", "Dominant", "Full", "Main", "Add-Diet", "Dom-Diet", "Full-Diet", "Diet")
    mainModels <- c("Additive", "Dominant", "Full", "Main")
    dietModels <- c("Add-Diet", "Dom-Diet", "Full-Diet", "Diet")
    dataType <- paste0(metab.names[i])
    
    dOut <- paste0(theDir, "plots/metab", i, "_", dataType, "/")
    if(!(file.exists(dOut))) dir.create(dOut)
    
    #All LOD values returned are compared to null (for main effects) and compared to full model (for diet effects) - thus get the incremental ones
    #All p-values shown are based on the incremental comparison
    results$LR2...Additive.incr <- results$LR...Additive
    results$LR2...Dominant.incr <- results$LR...Dominant - results$LR...Additive
    results$LR2...Full.incr <- results$LR...Full - results$LR...Dominant
    results$LR2...Main <- results$LR...Full
    results$LR2...D.Additive.incr <- results$LR...D.Additive - results$LR...Full
    results$LR2...D.Dominant.incr <- results$LR...D.Dominant - results$LR...D.Additive
    results$LR2...D.Full.incr <- results$LR...D.Full - results$LR...D.Dominant
    results$LR2...Diet <- results$LR...D.Full - results$LR...Full
    diffLOD <- names(results)[grepl("LR2", names(results))]
    
    #Peaks for the specific metabolite!
    mainCI.info <- paste0(dOut, '/mainCI_info_', dataType, '.csv')
    #write(paste0('Model Type,Metab Ind,Metab Name,Peak Ind,Peak Chr,Peak Pos (Mb),Lower Bound (Mb),Upper Bound (Mb),Length (Mb),p-value'), 
    #      file=mainCI.info, append=FALSE)
    dietCI.info <- paste0(dOut, '/dietCI_info_', dataType, '.csv')
    #write(paste0('Model Type,Metab Ind,Metab Name,Peak Ind,Peak Chr,Peak Pos (Mb),Lower Bound (Mb),Upper Bound (Mb),Length (Mb),p-value'), 
    #      file=dietCI.info, append=FALSE)
    for (n in 1:length(diffPVals)) {
      pVals.vec <- results[diffPVals[n]][[1]]
      
      #Find Peaks
      thres=-log10(0.05/length(pVals.vec))
      model.peaks <- pk_find(pVals.vec, thres=thres, model=fileNames[n])
      peakPos <- poslist[model.peaks,c("chr", "Ppos")]
      model.peaks.info <- cbind(rep(i,length(model.peaks)), rep(fileNames[n],length(model.peaks)), model.peaks, poslist[model.peaks,1], poslist[model.peaks,2]/1e6, 
                                results[model.peaks, diffPVals])
      colnames(model.peaks.info) <- c("Metab Ind", "Metab Name", "Peak Ind", "Peak Chr", "Peak Pos (Mb)", "P-Value - Add", "P-Value - Dom", "P-Value - Full", 
                                      "P-Value - Main", "P-Value - Add-D", "P-Value - Dom-D", "P-Value - Full-D", "P-Value - Diet")
      write(paste0(metab.names[i], '\t', diffPVals[n], '\t', length(model.peaks), '\t', paste(model.peaks, collapse = '\t')), file=metabPeaks.file, append=T)
      if (writePeaks.ind) {
        outFile <- paste0(dOut, "/peak_metab", i, "_", dataType, '_', fileNames[n], "_detailed.csv")
        write.table(model.peaks.info, file=outFile, sep=',')
      }
      
      if (findCIs) {
        #Find CI for each Peak
        qtldat <- data.frame(chr=CHR, Ppos=BP, Gpos=Gpos, LOD=results[diffLOD[n]][[1]])
        print("Finding CI!!")
        if (length(model.peaks) > 0) {
          #allCI <- matrix(0, ncol=6, nrow=nrow(peakPos))
          for (aPeak in 1:nrow(peakPos)) {
            #QTL Peak Info
            peak.qtl.ind <- model.peaks[aPeak]
            peak.qtl.chr <- peakPos[aPeak,1]
            peak.qtl.pos <- peakPos[aPeak,2]/1e6
            peak.pVal <- pVals.vec[model.peaks[aPeak]]
            
            #95% CI
            model.CI <- findCI(peakPos[aPeak,1], peakPos[aPeak,2], qtldat, method='LODdrop')
            CI.lower <- model.CI[1,2]/1e6
            CI.upper <- model.CI[2,2]/1e6
            CI.length <- CI.upper - CI.lower
            
            #Each Metabolite Peak File
            if (writePeaks.ind) {
              if (fileNames[n] %in% mainModels) {
                write(paste0(fileNames[n], ',', i, ',', dataType, ',', peak.qtl.ind, ',', peak.qtl.chr, ',', peak.qtl.pos, ',', CI.lower, ',', CI.upper, ',', CI.length, ',', peak.pVal), 
                      file=mainCI.info, append=TRUE)
              } else if (fileNames[n] %in% dietModels) {
                write(paste0(fileNames[n], ',', i, ',', dataType, ',', peak.qtl.ind, ',', peak.qtl.chr, ',', peak.qtl.pos, ',', CI.lower, ',', CI.upper, ',', CI.length, ',', peak.pVal), 
                      file=dietCI.info, append=TRUE)
              }
            }
            
            #All Metabolites Peak File
            #if (fileNames[n] %in% mainModels) {
            #  write(paste0(fileNames[n], ',', i, ',', dataType, ',', peak.qtl.ind, ',', peak.qtl.chr, ',', peak.qtl.pos, ',', CI.lower, ',', CI.upper, ',', CI.length, ',', peak.pVal), 
            #        file=allMainPeaks.file, append=TRUE)
            #} else if (fileNames[n] %in% dietModels) {
            #  write(paste0(fileNames[n], ',', i, ',', dataType, ',', peak.qtl.ind, ',', peak.qtl.chr, ',', peak.qtl.pos, ',', CI.lower, ',', CI.upper, ',', CI.length, ',', peak.pVal), 
            #        file=allDietPeaks.file, append=TRUE)
            #}
            
          }
          
          #Add p-value information to the appropriate matrix
          assign(paste0(diffPVals[n], '.matrix'), cbind(get(paste0(diffPVals[n], '.matrix')), c(dataType, pVals.vec)))
          
        } else {
          write.table("No Peaks", file=paste0(dOut, "/CI_metab", i, "_", fileNames[n], "2.csv"), row.names=F, sep=',')
        }
      }
      
      if (eachManhattan) {
        #Create Manhattan Plot
        print('Creating Manhattan Plots!!')
        if (length(model.peaks) > 0) {
          manhattanPlot(CHR, BP, pVals.vec, fileNames[n], dataType, dOut, peakPos=peakPos, thres=thres)
        } else {
          manhattanPlot(CHR, BP, pVals.vec, fileNames[n], dataType, dOut, peakPos=NA, thres=thres)
        }
      }
      
    }
  } else {
    print(paste0("Metabolite ", i, " does not have a file!!"))
  }
}

#save(P.Val...Additivve.matrix, file=paste0(theDir, "plots/add.matrix.rda"))
#save(P.Val...Dominant.matrix, file=paste0(theDir, 'plots/dom.matrix.rda'))
#save(P.Val...Full.matrix, file=paste0(theDir, 'plots/full.matrix.rda'))
#save(P.Val...Main.matrix, file=paste0(theDir, 'plots/main.matrix.rda'))
#save(P.Val...Add.D.matrix, file=paste0(theDir, 'plots/addDiet.matrix.rda'))
#save(P.Val...Dom.D.matrix, file=paste0(theDir, 'plots/domDiet.matrix.rda'))
#save(P.Val...Full.D.matrix, file=paste0(theDir, 'plots/fullDiet.matrix.rda'))
#save(P.Val...Diet.matrix, file=paste0(theDir, 'plots/diet.matrix.rda'))
