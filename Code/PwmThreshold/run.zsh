#!/bin/zsh

# Parameters
inLoc="/Users/egg/Projects/Gusmao_DeNovo/Results/Wagner_Filter/"
outLoc="/Users/egg/Projects/Gusmao_DeNovo/Results/Wagner_Filter/"
repList=( "PFM" )

# Execution
for rep in $repList
do
  python precomputeMotifThreshold.py $inLoc$rep"/" $outLoc$rep".fpr"
done


