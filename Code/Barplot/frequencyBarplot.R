
##########################################
### Initialization
##########################################

# Import
library(ggplot2)
library(reshape2)

# Input
args <- commandArgs(trailingOnly = TRUE)
graphWidth = as.numeric(args[1])
graphHeight = as.numeric(args[2])
titleSize = as.numeric(args[3])
tickSize = as.numeric(args[4])
frequencyTableFileName = args[5]
outputFileName = args[6]

# Barplot Parameters
barStat = "identity"
barPosition = "identity"
barAlpha = 0.7
scaleLabel = "Frequency"
scaleLabel1 = "Foreground"
scaleLabel2 = "Background"
xLab = "Transcription Factors"
yLab = "Frequency of Regions with Motif (%)"
barTitle = "Enriched Trancription Factors' Regions' Frequency"
deviceName = "pdf"
deviceDPI = 300


##########################################
### Heatmap Functions
##########################################

# Barplot
createBarplot <- function(dataFrame, outputFile){

  # Barplot
  datFMelt <- melt(dataFrame,id.vars = "FACTOR")
  bar = ggplot(data = datFMelt, aes(x = reorder(FACTOR, value), y = value, fill = variable))
  bar = bar + geom_bar(stat = barStat, position = barPosition, alpha = barAlpha)
  bar = bar + scale_fill_discrete(scaleLabel, labels=c(scaleLabel1, scaleLabel2))
  bar = bar + coord_flip()
  bar = bar + xlab(xLab) 
  bar = bar + ylab(yLab)  
  bar = bar + ggtitle(barTitle) 
  bar = bar + theme_bw()
  bar = bar + theme(plot.title = element_text(size=titleSize), axis.text.y = element_text(size=tickSize))

  # Saving graph
  ggsave(outputFile, plot=bar, device = deviceName, dpi = deviceDPI, width = graphWidth, height = graphHeight)

}

##########################################
### Execution
##########################################

# Reading Data Table
data = read.table(frequencyTableFileName, sep="\t", header=T, row.names=1)

# Barplot
outputFileNamePlot = paste(outputFileName, ".pdf", sep="")
createBarplot(data, outputFileNamePlot)



