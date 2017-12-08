
##########################################
### Initialization
##########################################

# Import
library(ggplot2)

# Barplot Parameters
graphWidth = 5
graphHeight = 5
deviceName = "pdf"
deviceDPI = 300
hmgb2Exp = 418.04

##########################################
### Heatmap Functions
##########################################

# Barplot
createDistribution <- function(dataFrame, intercPer, intercHmgb, outputFile){

  #bar = ggplot(data = dataFrame, aes(x = factor(dataFrame$SYMBOL, levels = dataFrame$SYMBOL[order(dataFrame$FC)]), y = dataFrame$FC, fill = dataFrame$ADJPVAL))
  #bar = bar + geom_bar(stat = barStat)
  #bar = bar + scale_fill_gradient2(low=lowColor, high=highColor, space=barSpace, limits=colorLimits, guide = colorGuide, name = colorName)
  #bar = bar + coord_flip()
  #bar = bar + xlab(xLab) 
  #bar = bar + ylab(yLab)  
  #bar = bar + ggtitle(barTitle) 
  #bar = bar + theme_bw()
  #bar = bar + geom_hline(yintercept=2, linetype = "dashed", colour = "green")
  #bar = bar + geom_hline(yintercept=0.5, linetype = "dashed", colour = "blue")
  #bar = bar + theme(plot.title = element_text(size=titleSize), axis.text.y = element_text(size=tickSize))

  # Barplot
  dist = ggplot(dataFrame, aes(x=exp))
  dist = dist + geom_histogram(aes(y=..density..), binwidth=.1, colour="black", fill="#54A7AA", alpha=0.8)
  dist = dist + geom_density(alpha=0.5, fill="#F47216")
  dist = dist + geom_vline(xintercept=intercPer, linetype = "dashed", colour = "#008000")
  dist = dist + geom_vline(xintercept=intercHmgb, linetype = "dashed", colour = "red")
  dist = dist + theme_bw()

  # Saving graph
  ggsave(outputFile, plot=dist, device = deviceName, dpi = deviceDPI, width = graphWidth, height = graphHeight)

}

##########################################
### Execution
##########################################

# Reading data frame
expressionTableFileName = "/Users/egg/Projects/Wagner_Methylation/Data/RNA_genes/RNA_genes_full_data.txt"
data = read.table(expressionTableFileName, sep="\t", header=T, na.strings=c("", "NA"))

# Subset data frame
newVec = log10(as.numeric(data[,12]))
quantExp = quantile(newVec, probs = c(0.01))[[1]]
newDataFr = data.frame(cond = factor(rep(c("WT"), each = length(newVec))), exp = newVec)

# Barplot
outputFileNamePlot = "/Users/egg/Projects/Wagner_Methylation/Results/ExpressionDistribution/RNA_genes_full_data.pdf"
createDistribution(newDataFr, quantExp, log10(hmgb2Exp), outputFileNamePlot)


