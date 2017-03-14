library(ggplot2)
theDir <- 'Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots'
setwd(theDir)
graph = T

metabScatter <- function(theMatrix) {
  sigIndices <- which(theMatrix$p.value>6)
  xCoord <- as.numeric(theMatrix$Ppos[sigIndices]/1e6)
  yCoord <- as.numeric(theMatrix$Metabolite.Index[sigIndices])
  chr <- as.character(theMatrix$chr[sigIndices])
  p.val <- theMatrix$p.value[sigIndices]
  groups <- sapply(p.val, function(x) if (x>10) {2} else if (x>8) {1} else {0})
  #groups <- yCoord
  theCoordinates <- cbind(xCoord, yCoord, chr, groups)
  return(theCoordinates)
}

paperFigures <- TRUE
allMain.sign <- read.csv(file='allMainPeaks.csv')
allDiet.sign <- read.csv(file='allDietPeaks.csv')

if (graph) {
  theModels <- c('Additive', 'Dominant', 'Full', 'Main', 'Add-Diet', 'Dom-Diet', 'Full-Diet', 'Diet')
  for (m in 1:length(theModels)) {
    if (m %in% c(1,2,3,4)) {
      model.data <- allMain.sign[which(allMain.sign$Model.Type==theModels[m]),]
    } else {
      model.data <- allDiet.sign[which(allDiet.sign$Model.Type==theModels[m]),]
    }
    
    model.groups <- sapply(model.data$p.value, function(x) if (x>10) {'10+'} else if (x>=8) {'8-10'} else {'5-8'})
    theColors <- c('blue', 'red', 'green', 'orange', 'black')
    color.categories <- c()
    for (x in 1:422) {
      remain <- x %% length(theColors)
      if (remain == 0) {
        remain = 5
      }
      color.categories <- c(color.categories, theColors[remain])
    }
    color.groups <- sapply(model.data$Metab.Ind, function(x) color.categories[x])
    theData <- data.frame(xAxis=as.double(model.data$QTL.Pos.Mb.), yAxis=as.double(model.data$Metab.Ind), chr=model.data$QTL.chr, grps=factor(model.groups, levels=c("5-8", "8-10", "10+")), colorGroups=color.groups)
    
    if (! paperFigures) {
      figureTitle <- paste0(theModels[m], " Model")
    } else {
      figureTitle <- NULL
    }
    
    ggplot(theData) + geom_point(aes(x=xAxis,y=yAxis,group=chr,color=grps,size=1.5)) + 
      facet_grid(.~chr,scales="fixed") + 
      labs(title=figureTitle,x="QTL Position (Mb)",y="Metabolite") + 
      scale_colour_manual(name="Significance", labels=c('5-8', '8-10', '10+'), breaks=c('5-8', '8-10', '10+'), values=c('blue', 'orange', 'red')) +
      theme(title=element_text(size=34), axis.title=element_text(size=24), axis.text.x=element_text(size=28, angle=45,hjust=0.5, vjust=0.5), axis.text.y=element_text(size=30), strip.text.x=element_text(size=30), 
            legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.height=unit(0.5, 'in')) + scale_x_continuous(breaks=c(0,10,20)) + guides(size=FALSE, colour=guide_legend(override.aes=list(size=4)))
    ggsave(filename=paste0("all_", theModels[m], "_sign_mB_paper_", paperFigures, ".pdf"), width=18, units='in')
    
    ggplot(theData) + geom_point(aes(x=xAxis,y=yAxis,group=chr,colour=as.character(colorGroups),size=1.5)) + 
      facet_grid(.~chr,scales="fixed") + 
      labs(title=figureTitle,x="QTL Position (Mb)",y="Metabolite") + 
      theme(title=element_text(size=34), axis.title=element_text(size=24), axis.text.x=element_text(size=28, angle=45,hjust=0.5, vjust=0.5), axis.text.y=element_text(size=30), strip.text.x=element_text(size=30), 
            legend.position='non') + scale_x_continuous(breaks=c(0,10,20)) + scale_colour_manual(values=c('blue', 'red', 'orange', 'darkgreen', 'purple'))
    ggsave(filename=paste0('all_', theModels[m], '_sign_mB_paper_', paperFigures, '_metabColor.pdf'), width=18, units='in')
  }
  
  merged.sign <- rbind(allMain.sign[which(allMain.sign$Model.Type == 'Main'),], allDiet.sign[which(allDiet.sign$Model.Type == 'Diet'),])
  #write.csv(merged.sign, file="genetic_gxe.csv", row.names=F, quote=F)
  theData <- data.frame(xAxis=merged.sign$QTL.Pos.Mb., yAxis=merged.sign$Metab.Ind, chr=merged.sign$QTL.chr, grps=merged.sign$Model.Type)
  
  if (! paperFigures) {
    figureTitle <- 'Genetic vs GxE'
  } else {
    figureTitle <- NULL
  }
  ggplot(theData) + geom_point(aes(x=xAxis,y=yAxis,group=chr,color=grps)) + 
    facet_grid(.~chr,scales='fixed') + 
    labs(title=figureTitle,x='QTL Position (Mb)',y='Metabolite') + 
    scale_colour_manual(name='Model Type', values=c('red', 'blue')) + 
    theme(title=element_text(size=34), axis.title=element_text(size=34), axis.text.x=element_text(size=30, angle=45,hjust=0.5,vjust=0.5), axis.text.y=element_text(size=30), strip.text.x=element_text(size=30), 
          legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.height=unit(0.5, 'in')) + scale_x_continuous(breaks=c(0,10,20))
  ggsave(filename=paste0('Genetic vs GxE_paper_', paperFigures, '.pdf'), width=18, units='in')
}
