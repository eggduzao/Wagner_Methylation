
##########################################
### Initialization
##########################################

# Import
library(ggplot2)

# Input
args <- commandArgs(trailingOnly = TRUE)
graphWidth = as.numeric(args[1])
graphHeight = as.numeric(args[2])
titleSize = as.numeric(args[3])
tickSize = as.numeric(args[4])
expressionTableFileName = args[5]
geneSymbolsTable = args[6]
outputFileName = args[7]

# Barplot Parameters
barStat = "identity"
lowColor = "white"
highColor = "red"
barSpace = "Lab"
minLimCol = 0
maxLimCol = 3
colorLimits = c(minLimCol, maxLimCol)
colorGuide = "colourbar"
colorName = "P-value\n(-log10)"
xLab = "Transcription Factors"
yLab = "Gene Expression Fold Change"
barTitle = "Enriched Trancription Factors' Gene Expression"
deviceName = "pdf"
deviceDPI = 300

##########################################
### Heatmap Functions
##########################################

# Barplot
createBarplot <- function(dataFrame, outputFile){

  # Barplot
  bar = ggplot(data = dataFrame, aes(x = factor(dataFrame$SYMBOL, levels = dataFrame$SYMBOL[order(dataFrame$FC)]), y = dataFrame$FC, fill = dataFrame$ADJPVAL))
  bar = bar + geom_bar(stat = barStat)
  bar = bar + scale_fill_gradient2(low=lowColor, high=highColor, space=barSpace, limits=colorLimits, guide = colorGuide, name = colorName)
  bar = bar + coord_flip()
  bar = bar + xlab(xLab) 
  bar = bar + ylab(yLab)  
  bar = bar + ggtitle(barTitle) 
  bar = bar + theme_bw()
  bar = bar + geom_hline(yintercept=2, linetype = "dashed", colour = "green")
  bar = bar + geom_hline(yintercept=0.5, linetype = "dashed", colour = "blue")
  bar = bar + theme(plot.title = element_text(size=titleSize), axis.text.y = element_text(size=tickSize))

  # Saving graph
  ggsave(outputFile, plot=bar, device = deviceName, dpi = deviceDPI, width = graphWidth, height = graphHeight)

}

##########################################
### Execution
##########################################

# Reading gene symbols
geneSymbols = read.table(geneSymbolsTable, sep="\t", header=F)
geneSymbols = as.character(geneSymbols[,1])

# Reading Data Table
data = read.table(expressionTableFileName, sep="\t", header=T, row.names=2, na.strings=c("", "NA"))
data = data[geneSymbols,c(1,3,4)]
data[,1] = rownames(data)
colnames(data) = c("SYMBOL", "ADJPVAL", "FC")
data = data[complete.cases(data), ]

# P-value in -log10 scale
data[,2] = -log10(data[,2])

# Subset p-value to be not greater than the maximum color plot
data[data[,2]>maxLimCol,2] = maxLimCol

# Barplot
outputFileNamePlot = paste(outputFileName, ".pdf", sep="")
createBarplot(data, outputFileNamePlot)



