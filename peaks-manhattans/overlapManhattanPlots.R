library(DSPRqtl)
library(ggplot2)
library(hash)
library(RColorBrewer)
library(reshape2)

theDir <- 'Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/'
setwd(theDir)
mainSubset <- TRUE
dietSubset <- TRUE
paperFigures <- TRUE

manhattanPlot_overlay <- function(chr.vec, bp.vec, pVal.matrix, modelType, thres=6) {
  removeNA <- which(is.na(pVal.matrix), arr.ind=T)[,1]
  if (length(removeNA) > 0) {
    pVal.matrix <- pVal.matrix[-removeNA,]
    chr.vec <- chr.vec[-removeNA]
    bp.vec <- bp.vec[-removeNA]
  }
  
  theData <- data.frame(bp=bp.vec, chr=chr.vec, pVals=apply(pVal.matrix, 2, as.numeric))
  colnames(theData) <- c('bp', 'chr', colnames(pVal.matrix))
  theData.melt <- melt(theData, id.vars=1:2, measure.vars=3:ncol(theData))
  
  if (! paperFigures) {
    if (modelType == 'add') {
      figureTitle <- 'Genetic Additive Effect mQTLs'
    } else if (modelType == 'dom') {
      figureTitle <- 'Genetic Dominance Effect mQTLs'
    } else if (modelType == 'full') {
      figureTitle <- 'Genetic Full Effect mQTLs'
    } else if (modelType == 'main') {
      figureTitle <- 'Genetic Main Effect mQTLs'
    } else if (modelType == 'addDiet') {
      figureTitle <- 'Add x Diet Effect mQTLs'
    } else if (modelType == 'domDiet') {
      figureTitle <- 'Dom x Diet Effect mQTLs'
    } else if (modelType == 'fullDiet') {
      figureTitle <- 'Full x Diet Effect mQTLs'
    } else if (modelType == 'diet') {
      figureTitle <- 'Genetic x Diet Effect mQTLs'
    }
  } else {
    figureTitle <- NULL
  }
  
  y.limit <- ceiling(max(theData.melt$value))
  thePlot <- ggplot(theData.melt) + geom_point(aes(x=bp,y=value,group=chr, colour=variable)) + 
    facet_grid(.~chr,scales="fixed") + geom_hline(yintercept = thres, col='black', size=1.5) +
    labs(title=figureTitle,x="QTL Position (Mb)",y="-log10(p)") + scale_y_continuous(limits = c(5, y.limit)) + 
    guides(colour=FALSE) + theme(title=element_text(size=34), axis.title=element_text(size=34, vjust=1), axis.text.x=element_text(size=30, angle=45, hjust=0.5, vjust=0.5), axis.text.y=element_text(size=30), strip.text.x=element_text(size=30))
  ggsave(filename=paste0(modelType, '.plot.paper_', paperFigures, '.pdf'), plot=thePlot, width=10, units='in')
  
  if (mainSubset & modelType=='main') {
    metaboliteSubsets <- c('citrulline', 'glycolic_acid', 'lactic_acid', 'N_epsilon_trimethyllysine', 'palmitoleic_acid', 
                           'pantothenic_acid', 'shikimic_acid')
    theData.subset <- theData[,c(1, 2, which(colnames(theData) %in% metaboliteSubsets))]
    theData.subset.melt <- melt(theData.subset, id.vars=1:2, measure.vars=3:ncol(theData.subset))
    y.limit <- ceiling(max(theData.subset.melt$value))
    
    thePlot <- ggplot(theData.subset.melt) + geom_point(aes(x=bp,y=value,group=chr, colour=variable)) + 
      facet_grid(.~chr,scales="fixed") + geom_hline(yintercept = thres, col='black', size=1.25) +
      labs(title=figureTitle,x="QTL Position (Mb)",y="-log10(p)") + scale_y_continuous(limits = c(5, y.limit)) + 
      theme(title=element_text(size=24), axis.title=element_text(size=24), axis.text.x=element_text(size=20, angle=45, hjust=0.5, vjust=0.5), axis.text.y=element_text(size=20), strip.text.x=element_text(size=20), 
            legend.text=element_text(size=18), legend.title=element_text(size=20))
    thePlot + scale_fill_discrete(name='Metabolite', breaks=metaboliteSubsets, label=metaboliteSubsets)
    ggsave(filename=paste0(modelType, '.subset.plot.paper_', paperFigures, '.pdf'), plot=thePlot, width=10, units='in')
  }
  
  if (dietSubset & modelType=='diet') {
    metaboliteSubsets <- c('alanine', 'cytidine_5_monophosphate_NIST', 'hydroxycarbamate_NIST', 'inosine', 'lactic_acid', 'myo_inositol', 
                           'urea')
    theData.subset <- theData[,c(1, 2, which(colnames(theData) %in% metaboliteSubsets))]
    theData.subset.melt <- melt(theData.subset, id.vars=1:2, measure.vars=3:ncol(theData.subset))
    y.limit <- ceiling(max(theData.subset.melt$value))
    
    thePlot <- ggplot(theData.subset.melt) + geom_point(aes(x=bp,y=value,group=chr, colour=variable)) + 
      facet_grid(.~chr,scales="free") + geom_hline(yintercept = thres, col='black', size=1.25) +
      labs(title=figureTitle,x="QTL Position (Mb)",y="-log10(p)") + scale_y_continuous(limits = c(5, y.limit)) + 
      theme(title=element_text(size=24), axis.title=element_text(size=24), axis.text.x=element_text(size=20, angle=45, hjust=0.5, vjust=0.5), axis.text.y=element_text(size=20), strip.text.x=element_text(size=20), 
            legend.text=element_text(size=18), legend.title=element_text(size=20))
    thePlot + scale_fill_discrete(name='Metabolite', breaks=metaboliteSubsets, label=metaboliteSubsets)
    ggsave(filename=paste0(modelType, '.subset.plot.paper_', paperFigures, '.pdf'), plot=thePlot, width=10, units='in')
  }
  
}

models.fileNames <- c('add', 'dom', 'full', 'main', 'addDiet', 'domDiet', 'fullDiet', 'diet')
models.objNames <- c("P.Val...Additivve.matrix", "P.Val...Dominant.matrix", "P.Val...Full.matrix", "P.Val...Main.matrix", "P.Val...Add.D.matrix", "P.Val...Dom.D.matrix", 
                     "P.Val...Full.D.matrix", "P.Val...Diet.matrix")
for (m in 1:length(models.fileNames)) {
  load(paste0(models.fileNames[m], '.matrix.rda'))
  theMatrix <- get(models.objNames[m])
  colnames(theMatrix) <- theMatrix[1,]
  theMatrix <- theMatrix[-1,]
  
  chr.vec <- poslist[,1]
  bp.vec <- poslist[,2]/1e6
  thres <- -log10(0.05/nrow(theMatrix))
  thePlot <- manhattanPlot_overlay(chr.vec, bp.vec, theMatrix, models.fileNames[m], thres=thres)
}
