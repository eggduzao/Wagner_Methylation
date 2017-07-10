
# Import
import os
import sys

# Input
expFileName = "/home/egg/Projects/Wagner_Methylation/Data/RNA_genes/RNA_genes_full_data.txt"
aliasFileName = "/home/egg/Projects/Wagner_Methylation/Data/genome/alias_human.txt"
genesFileName = "/home/egg/Projects/Wagner_Methylation/Data/genome/genes_RefSeq_hg19.bed"
dnasePeakFileName = "/home/egg/Projects/Gusmao_DeNovo/Wagner_Data/MSC_OpenChromatin/peaks/DNase_MSC_filtered_peaks.bed"
outputFileNameUpReg = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_up.bed"
outputFileNameDownReg = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_down.bed"

# Parameters
geneBodyAdd = 20
promoterSize = 2000
pValueCutoff = 0.05
upFCCutoff = 2
downFCCutoff = 0.5

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

# Creating gene promoter positions dictionary
# GENE -> [[chrom, pos1, pos2], ...]
genesDict = dict()
genesFile = open(genesFileName,"r")
for line in genesFile:
  ll = line.strip().split("\t")
  chrom = ll[0]
  geneN = ll[3].upper()
  if(ll[5] == "+"):
    p1 = str(max(0,int(ll[1])-promoterSize))
    p2 = str(int(ll[1])+geneBodyAdd)
  else:
    p1 = str(int(ll[2])-geneBodyAdd)
    p2 = str(int(ll[2])+promoterSize)
  try: genesDict[aliasDict[geneN]].append([chrom,p1,p2])
  except Exception:
    try: genesDict[aliasDict[geneN]] = [[chrom,p1,p2]]
    except Exception: genesDict[geneN] = [[chrom,p1,p2]]

# Iterating on methylation table
toRemove = []
expFile = open(expFileName,"r")
for e in range(0,1): expFile.readline()
outputFileNameUpRegTemp = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_up_temp.bed"
outputFileNameDownRegTemp = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_down_temp.bed"
toRemove.append(outputFileNameUpRegTemp); toRemove.append(outputFileNameDownRegTemp)
outputFileUpRegTemp = open(outputFileNameUpRegTemp,"w")
outputFileDownRegTemp = open(outputFileNameDownRegTemp,"w")
for line in expFile:

  # Initialization
  ll = line.strip().split("\t")
  try: gene = aliasDict[ll[1].upper()]
  except Exception: methGene = "NOGENE"
  try: assocGeneList = genesDict[gene]
  except Exception: assocGeneList = []
  try: pvalue = float(ll[3])
  except Exception: continue
  try: foldchange = float(ll[4])
  except Exception: continue

  # Filter1
  if(pvalue >= pValueCutoff): continue

  # Filter2
  if(foldchange >= upFCCutoff): # Upregulated
    for vec in assocGeneList: outputFileUpRegTemp.write("\t".join([vec[0],vec[1],vec[2]])+"\n")
  elif(foldchange <= downFCCutoff): # Downregulated
    for vec in assocGeneList: outputFileDownRegTemp.write("\t".join([vec[0],vec[1],vec[2]])+"\n")

expFile.close()
outputFileUpRegTemp.close()
outputFileDownRegTemp.close()

# Sorting output files
outputFileNameUpRegTempSort = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_up_sort.bed"
outputFileNameDownRegTempSort = "/home/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/genes_down_sort.bed"
toRemove.append(outputFileNameUpRegTempSort); toRemove.append(outputFileNameDownRegTempSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameUpRegTemp+" | uniq > "+outputFileNameUpRegTempSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameDownRegTemp+" | uniq > "+outputFileNameDownRegTempSort)

# Intersecting with open chromatin region
os.system("intersectBed -a "+outputFileNameUpRegTempSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameUpReg)
os.system("intersectBed -a "+outputFileNameDownRegTempSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameDownReg)

# Removing files
for e in toRemove: os.system("rm "+e)


