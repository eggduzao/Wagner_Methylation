
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)

# Input
args <- commandArgs(trailingOnly = TRUE)
parCexMain = as.numeric(args[1])
graphWidth = as.numeric(args[2])
graphHeight = as.numeric(args[3])
heatmapCexCol = as.numeric(args[4])
heatmapCexRow = as.numeric(args[5])
heatmapCexRow2 = as.numeric(args[6])
heatmapLhei1 = as.numeric(args[7])
heatmapLhei2 = as.numeric(args[8])
graphWidthD = as.numeric(args[9])
graphHeightD = as.numeric(args[10])
heatmapCexColD = as.numeric(args[11])
heatmapCexRowD = as.numeric(args[12])
heatmapCexRow2D = as.numeric(args[13])
heatmapLhei1D = as.numeric(args[14])
heatmapLhei2D = as.numeric(args[15])
title1 = args[16]
title1D = args[17]
title2 = args[18]
title2D = args[19]
inputTableFileName = args[20]
outputFileName = args[21]

# General Parameters
pvalue1 = 0.05
pvalue2 = 0.05
ltpvalue1 = -log10(pvalue1) - 0.05
title1 = gsub("_", " ", title1)
title1D = gsub("_", " ", title1D)
title2 = gsub("_", " ", title2)
title2D = gsub("_", " ", title2D)
title1 = gsub("NN", "\n", title1)
title1D = gsub("NN", "\n", title1D)
title2 = gsub("NN", "\n", title2)
title2D = gsub("NN", "\n", title2D)
title1 = gsub("PP", "(", title1)
title1D = gsub("PP", "(", title1D)
title2 = gsub("PP", "(", title2)
title2D = gsub("PP", "(", title2D)
title1 = gsub("QQ", ")", title1)
title1D = gsub("QQ", ")", title1D)
title2 = gsub("QQ", ")", title2)
title2D = gsub("QQ", ")", title2D)

# Regular Heatmap Parameters
clusteringDistance = "euclidean"
clusteringMethod = "ward"
colorScheme = colorRampPalette(brewer.pal(9, 'Reds'))(149)
hmbreaks = c(seq(0,5,length=150))
colorScheme[hmbreaks[-length(hmbreaks)]<ltpvalue1] = "gray"
parCexAxis = 1.0
heatmapStrCol = 0
heatmapAdjCol = c(0.5,1)
heatmapMargins = c(5,10)
heatmapSepWidth = c(2,2)
heatmapLhei = c(heatmapLhei1, heatmapLhei2)
heatmapSepColor = "black"
heatmapKeySize = 2.5
heatmapTrace = "none"
heatmapTraceCol = NA
heatmapDensityInfo = "none"
heatmapRowv = TRUE
heatmapColv = FALSE

# Differential Heatmap Parameters
colorSchemeD = colorRampPalette(brewer.pal(9, 'RdBu'))(149)
hmbreaksD = c(seq(-5,5,length=150))
colorSchemeD[hmbreaksD[-length(hmbreaksD)] > -1 & 1 > hmbreaksD[-length(hmbreaksD)]] = "gray"
parCexAxisD = 1.0
heatmapLabColD = NA
heatmapMarginsD = c(5,10)
heatmapSepWidthD = c(2,2)
heatmapLheiD = c(heatmapLhei1D, heatmapLhei2D)
heatmapSepColorD = "black"
heatmapKeySizeD = 2.5
heatmapTraceD = "none"
heatmapTraceColD = NA
heatmapDensityInfoD = "none"
heatmapDendrogramD = "none"
heatmapRowvD = FALSE
heatmapColvD = FALSE

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap <- function(data, myTitle, heatmapCexRow, outputFile){

  # Clustering
  myDist = function(p1) dist(p1, method = clusteringDistance)
  myHclust = function(p2) hclust(p2, method = clusteringMethod)

  # Initializing graph
  postscript(outputFile, width = graphWidth, height = graphHeight, horizontal = FALSE, paper = "special")
  par(cex.axis = parCexAxis, cex.main = parCexMain)

  # Heatmap
  heatmap.2(data, col=colorScheme, breaks = hmbreaks, distfun = myDist, hclustfun = myHclust,
            main = myTitle, cexCol = heatmapCexCol, cexRow = heatmapCexRow, srtCol = heatmapStrCol, adjCol = heatmapAdjCol,
            margins = heatmapMargins, sepwidth = heatmapSepWidth, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, Rowv = heatmapRowv, Colv = heatmapColv)

  # Closing graph
  dev.off()

  # Remove eps
  system(paste("epstopdf \"",outputFile,"\"",sep=""))
  system(paste("rm \"",outputFile,"\"",sep=""))

}

# Differential Heatmap
createHeatmapDiff <- function(data, myTitle, heatmapCexRow, outputFile){

  # Initializing graph
  postscript(outputFile, width = graphWidthD, height = graphHeightD, horizontal = FALSE, paper = 'special')
  par(cex.axis = parCexAxisD, cex.main = parCexMain)

  # Heatmap
  heatmap.2(data, col = colorSchemeD, breaks = hmbreaksD,
            main = myTitle, cexCol = heatmapCexColD, cexRow = heatmapCexRow, labCol = heatmapLabColD, 
            margins = heatmapMarginsD, sepwidth = heatmapSepWidthD, lhei = heatmapLheiD, sepcolor = heatmapSepColorD, keysize = heatmapKeySizeD,
            trace = heatmapTraceD, tracecol = heatmapTraceColD, density.info = heatmapDensityInfoD,
            dendrogram = heatmapDendrogramD, Rowv = heatmapRowvD, Colv = heatmapColvD)

  # Closing graph
  dev.off()

  # Remove eps
  system(paste("epstopdf \"",outputFile,"\"",sep=""))
  system(paste("rm \"",outputFile,"\"",sep=""))

}

# Printing genes and motifs
printGenesAndMotifs <- function(mydata, outputFileGenes, outputFileMotifs){

  # Genes
  outTemp = "genestemp.txt"
  toWrite = c()
  for(r in rownames(mydata)){
    if(grepl("::", r)){
      ss = strsplit(r, "::")[[1]]
      for(s in ss){ toWrite = c(toWrite, s) }
    } else{ toWrite = c(toWrite, r) }
  }
  write(toWrite, file = outTemp, ncolumns = 1, append = FALSE, sep="")
  system(paste("sort ",outTemp," | uniq > ",outputFileGenes,sep=""))
  system(paste("rm ",outTemp,sep=""))

  # Motifs
  write(as.vector(mydata$MOTIF), file = outTemp, ncolumns = 1, append = FALSE, sep="")
  system(paste("sort ",outTemp," | uniq > ",outputFileMotifs,sep=""))
  system(paste("rm ",outTemp,sep=""))

}

##########################################
### Execution
##########################################

data = read.table(inputTableFileName, sep="\t", header=T, row.names=1)
dataN = data[,c(1,2)]
minV = min(dataN[dataN > 0]) # Calculating the minimum value greater than 0.0 

# Filter 1 = Keep only the lines in which there is at least one value < p-value1
dataFilt1 = data[apply(data[,c(1,2)], 1, function(x) sum(x < pvalue1)) > 0, ]

# Filter 2 = Keep only the lines in which there is at least one value > p-value2
dataFilt2 = dataFilt1[apply(dataFilt1[,c(1,2)], 1, function(x) sum(x > pvalue2)) > 0, ]

# Print genes and motifs
outputFileName1G = paste(outputFileName,"_1_genes.txt",sep="")
outputFileName1M = paste(outputFileName,"_1_motifs.txt",sep="")
printGenesAndMotifs(dataFilt1,outputFileName1G,outputFileName1M)
outputFileName2G = paste(outputFileName,"_2_genes.txt",sep="")
outputFileName2M = paste(outputFileName,"_2_motifs.txt",sep="")
printGenesAndMotifs(dataFilt2,outputFileName2G,outputFileName2M)

# Converting dataFilts to matrices
dataFilt1 = as.matrix(dataFilt1[,c(1,2)])
dataFilt2 = as.matrix(dataFilt2[,c(1,2)])

# All values X are now -log10(X+pseudocount), where pseudocount is the minimum value > 0.
dataFilt1 = apply(dataFilt1, 1:2, function(x) -log10(x+minV)) # Converting values
dataFilt2 = apply(dataFilt2, 1:2, function(x) -log10(x+minV)) # Converting values

# Creating Regular Heatmap
outputFileName1 = paste(outputFileName,"_1.eps",sep="")
createHeatmap(dataFilt1, title1, heatmapCexRow, outputFileName1)
outputFileName2 = paste(outputFileName,"_2.eps",sep="")
createHeatmap(dataFilt2, title2, heatmapCexRow2, outputFileName2)

# Calculating Differences
dataFilt1D = dataFilt1[,1,drop=FALSE] - dataFilt1[,2,drop=FALSE]
dataFilt2D = dataFilt2[,1,drop=FALSE] - dataFilt2[,2,drop=FALSE]

# Making diff double for heatmap
dataFilt1D = cbind(dataFilt1D,dataFilt1D)
dataFilt2D = cbind(dataFilt2D,dataFilt2D)

# Sorting diff matrices
dataFilt1D = dataFilt1D[order(dataFilt1D[,1]),]
dataFilt2D = dataFilt2D[order(dataFilt2D[,1]),]

# Creating Differential Heatmap
outputFileName1 = paste(outputFileName,"_1D.eps",sep="")
createHeatmapDiff(dataFilt1D, title1D, heatmapCexRowD, outputFileName1)
outputFileName2 = paste(outputFileName,"_2D.eps",sep="")
createHeatmapDiff(dataFilt2D, title2D, heatmapCexRow2D, outputFileName2)


