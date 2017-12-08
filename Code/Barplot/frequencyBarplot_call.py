
##########################################
### Initialization
##########################################

# Import
import os
import sys

# Input
graphWidth = sys.argv[1]
graphHeight = sys.argv[2]
titleSize = sys.argv[3]
tickSize = sys.argv[4]
frequencyFileNameUp = sys.argv[5]
frequencyFileNameDown = sys.argv[6]
geneSymbolsTable = sys.argv[7]
outputFileName = sys.argv[8]

# Initialization
toRemove = []
filePath = "/".join(os.path.realpath(__file__).split("/")[:-1])+"/"

# Function
def updateFrequencyDict(freqDict, geneSymbVec, frequencyFileName):
  frequencyFile = open(frequencyFileName,"r")
  for line in frequencyFile:
    ll = line.strip().split("\t")
    k = ll[0]
    if(k in geneSymbVec):
      pv = float(ll[1])
      cpv = float(ll[2])
      freq = ll[7]
      freq = freq[:-1]
      freqb = ll[8]
      freqb = freqb[:-1]
      if(k[:3] == "MA0" or k[:3] == "MA1"):
        factor = k.split(".")[2].upper()
      elif(k[:3] == "UP0"):
        factor = k.split("_")[2].upper()
      elif("filter" in k):
        factor = k.split("_")[2].upper()
      else:
        factor = k.split("_")[0].upper()
      try: 
        if(freqDict[k][4] > pv): freqDict[k] = [factor,factor,freq,freqb,pv]
      except Exception: freqDict[k] = [factor,factor,freq,freqb,pv]
  frequencyFile.close()

# Reading gene sumbols
geneSymbolsFile = open(geneSymbolsTable,"r")
geneSymbols = []
for line in geneSymbolsFile: geneSymbols.append(line.strip())
geneSymbolsFile.close()

# Fetching frequencyDict
frequencyDict = dict()
updateFrequencyDict(frequencyDict, geneSymbols, frequencyFileNameUp)
updateFrequencyDict(frequencyDict, geneSymbols, frequencyFileNameDown)

# Writing table
tempOutTable = "./freq.txt"; toRemove.append(tempOutTable)
tempOutFile = open(tempOutTable,"w")
tempOutFile.write("\t".join(["FACTOR", "FREQ", "BACKFREQ"])+"\n")
for k in sorted(frequencyDict.keys()): tempOutFile.write("\t".join(frequencyDict[k][:4])+"\n")
tempOutFile.close()

# Appying R
os.system("R CMD BATCH '--args '"+graphWidth+"' '"+graphHeight+"' '"+titleSize+"' '"+tickSize+"' '"+tempOutTable+"' '"+outputFileName+" "+filePath+"frequencyBarplot.R ./frequencyBarplot.Rout")
toRemove.append("./frequencyBarplot.Rout"); toRemove.append(".RData")

# Removing files
for e in toRemove: os.system("rm "+e)


