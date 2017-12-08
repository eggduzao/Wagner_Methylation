
#################################################
# Input
#################################################

# Import
import os
import sys

# Input
aliasFileName = sys.argv[1]
colName1 = sys.argv[2]
colName2 = sys.argv[3]
inputResultFileName1 = sys.argv[4]
inputResultFileName2 = sys.argv[5]
expressionFileName1 = sys.argv[6]
if(expressionFileName1 == "."): expressionFileName1 = None
expressionFileName2 = sys.argv[7]
if(expressionFileName2 == "."): expressionFileName2 = None
outputFileName = sys.argv[8]

# Parameters
pValueCutoff = 0.05
upFCCutoff = 2
downFCCutoff = 0.5

#################################################
# Fetch results of motif matching
#################################################

# Create factor dictionary
factorDict = dict()
inputResultFile1 = open(inputResultFileName1,"r")
inputResultFile1.readline()
for line in inputResultFile1:
  ll = line.strip().split("\t")
  factorDict[ll[0]] = [ll[1]]
inputResultFile1.close()
inputResultFile2 = open(inputResultFileName2,"r")
inputResultFile2.readline()
for line in inputResultFile2:
  ll = line.strip().split("\t")
  factorDict[ll[0]].append(ll[1])
inputResultFile2.close()

# Create a single-factor dictionary
sfactorDict = dict()
for k in factorDict.keys():
  if(k[:3] == "MA0" or k[:3] == "MA1"):
    factor = k.split(".")[2].upper()
  elif(k[:3] == "UP0"):
    factor = k.split("_")[2].upper()
  elif("filter" in k):
    factor = k.split("_")[2].upper()
  else:
    factor = k.split("_")[0].upper()
  n = factorDict[k][0]
  try:
    if(float(sfactorDict[factor][0]) > float(n)): sfactorDict[factor] = factorDict[k]+[k]
  except Exception: sfactorDict[factor] = factorDict[k]+[k]

#################################################
# Fetching gene expression
#################################################

# Creating alias dictionary
# GENE_SYMBOL -> GENE_ID
aliasDict = dict()
aliasFile = open(aliasFileName,"r")
for line in aliasFile:
  ll = line.strip().split("\t")
  ensGene = ll[0].upper()
  symGene = ll[1].upper()
  lll = ll[2].split("&")
  aliasDict[ensGene] = ensGene
  aliasDict[symGene] = ensGene
  for e in lll: aliasDict[e.upper()] = ensGene
aliasFile.close()

# Iterating on expression table
toRemove = []
try: expFile = open(expressionFileName1,"r")
except Exception:
  try: expFile = open(expressionFileName2,"r")
  except Exception: expFile = None
if(expFile):
  for e in range(0,1): expFile.readline()
  expDict = dict()
  for line in expFile:

    # Initialization
    ll = line.strip().split("\t")
    try: gene = aliasDict[ll[1].upper()]
    except Exception: gene = ll[1].upper()
    try: pvalue = float(ll[3])
    except Exception: continue
    try: foldchange = float(ll[4])
    except Exception: continue

    # Filter1
    if(pvalue >= pValueCutoff): continue

    # Filter2
    if(foldchange >= upFCCutoff or foldchange <= downFCCutoff):
      expDict[gene] = True

  expFile.close()

#################################################
# Output
#################################################

# Writing table
outputFile = open(outputFileName,"w")
outputFile.write("\t".join(["GENE", colName1, colName2, "MOTIF"])+"\n")
for k in sorted(sfactorDict.keys()):
  k1 = k.split("_F")[0]
  k2 = k1.split("::")
  kk = []
  for e in k2:
    kk.append(e)
    kk.append(e.split("-")[0])
  vec = [k]
  if(expressionFileName1):
    flag = False
    for e in kk:
      try: expDict[aliasDict[e]]
      except Exception:
        try: expDict[e]
        except Exception: continue
      flag = True
      vec.append(sfactorDict[k][0])
      break
    if(not flag): vec.append("1.0")
  else: vec.append(sfactorDict[k][0])
  if(expressionFileName2):
    flag = False
    for e in kk:
      try: expDict[aliasDict[e]]
      except Exception:
        try: expDict[e]
        except Exception: continue
      flag = True
      vec.append(sfactorDict[k][1])
      break
    if(not flag): vec.append("1.0")
  else: vec.append(sfactorDict[k][1])
  vec.append(sfactorDict[k][2])
  outputFile.write("\t".join(vec)+"\n")
outputFile.close()


