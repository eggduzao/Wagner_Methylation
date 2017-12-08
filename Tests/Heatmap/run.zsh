#!/usr/bin/env zsh

# Create tables
aliasFileName="/Users/egg/rgtdata/hg19/alias_human.txt"
colName1="Hypermethylated"
colName2="Hypomethylated"
inputResultFileName1="./ME1.txt"
inputResultFileName2="./ME2.txt"
expressionFileName1="./RNA.txt"
expressionFileName2="./RNA.txt"
outputFileName="./table.txt"
python /Users/egg/Projects/Wagner_Methylation/Code/Heatmap/createHeatmapTables.py $aliasFileName $colName1 $colName2 $inputResultFileName1 $inputResultFileName2 $expressionFileName1 $expressionFileName2 $outputFileName

# Create hetmaps
title1="Up-_and_Down-regulated_Genes_PPV1QQNN-log10PPp-valueQQ"
title1D="Title_1_(D)"
title2="Title_2"
title2D="Title_2_(D)"
inputTableFileName="./table.txt"
outputFileName="./heatmap"
R CMD BATCH '--args '$title1' '$title1D' '$title2' '$title2D' '$inputTableFileName' '$outputFileName /Users/egg/Projects/Wagner_Methylation/Code/Heatmap/heatmap.R ./heatmap.Rout


