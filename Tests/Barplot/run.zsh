#!/usr/bin/env zsh

# Create tables
expressionTableFileName="./RNA.txt"
geneSymbols="./genes.txt"
outputFileName="./expression"
R CMD BATCH '--args '$expressionTableFileName' '$geneSymbols' '$outputFileName /Users/egg/Projects/Wagner_Methylation/Code/Barplot/expressionBarplot.R ./expressionBarplot.Rout

# Create tables
frequencyFileName1="./ME_20_1.txt"
frequencyFileName2="./ME_20_2.txt"
geneSymbols="./motifs.txt"
outputFileName="./frequency"
python /Users/egg/Projects/Wagner_Methylation/Code/Barplot/frequencyBarplot_call.py $frequencyFileName1 $frequencyFileName2 $geneSymbols $outputFileName


