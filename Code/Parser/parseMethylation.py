
# Import
import os
import sys

# Input
methFileName = "/Users/egg/Projects/Wagner_Methylation/Data/DNAm/DNAm_full_data.txt"
aliasFileName = "/Users/egg/rgtdata/hg19/alias_human.txt"
genesFileName = "/Users/egg/rgtdata/hg19/genes_RefSeq_hg19.bed"
dnasePeakFileName = "/Users/egg/Projects/Gusmao_DeNovo/Data/Wagner_Data/MSC_OpenChromatin/peaks/DNase_MSC_filtered_peaks.bed"
outLoc = "/Users/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/"
outputFileNameHyperMeth = outLoc+"hyper_meth.bed"
outputFileNameHypoMeth = outLoc+"hypo_meth.bed"
outputFileNameHyperMethGenes = outLoc+"hyper_meth_genes.bed"
outputFileNameHypoMethGenes = outLoc+"hypo_meth_genes.bed"

# Parameters
geneBodyAdd = 20
promoterSize = 2000
pValueCutoff = 0.01
hyperCutoff = 0.05
hypoCutoff = -0.05
cpgExtend = 20

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
# GENE -> [[pos1, pos2], ...]
genesDict = dict()
genesFile = open(genesFileName,"r")
for line in genesFile:
  ll = line.strip().split("\t")
  geneN = ll[3].upper()
  if(ll[5] == "+"):
    p1 = str(max(0,int(ll[1])-promoterSize))
    p2 = str(int(ll[1])+geneBodyAdd)
  else:
    p1 = str(int(ll[2])-geneBodyAdd)
    p2 = str(int(ll[2])+promoterSize)
  try: genesDict[aliasDict[geneN]].append([p1,p2])
  except Exception:
    try: genesDict[aliasDict[geneN]] = [[p1,p2]]
    except Exception: genesDict[geneN] = [[p1,p2]]

# Iterating on methylation table
toRemove = []
methFile = open(methFileName,"r")
for e in range(0,2): methFile.readline()
outputFileNameHyperMethTemp = outLoc+"hyper_meth_temp.bed"
outputFileNameHypoMethTemp = outLoc+"hypo_meth_temp.bed"
outputFileNameHyperMethGenesTemp = outLoc+"hyper_meth_genes_temp.bed"
outputFileNameHypoMethGenesTemp = outLoc+"hypo_meth_genes_temp.bed"
toRemove.append(outputFileNameHyperMethTemp); toRemove.append(outputFileNameHypoMethTemp)
toRemove.append(outputFileNameHyperMethGenesTemp); toRemove.append(outputFileNameHypoMethGenesTemp)
outputFileHyperMeth = open(outputFileNameHyperMethTemp,"w")
outputFileHypoMeth = open(outputFileNameHypoMethTemp,"w")
outputFileHyperMethGenes = open(outputFileNameHyperMethGenesTemp,"w")
outputFileHypoMethGenes = open(outputFileNameHypoMethGenesTemp,"w")
for line in methFile:

  # Initialization
  ll = line.strip().split("\t")
  try: methGene = aliasDict[ll[1].upper()]
  except Exception: methGene = "NOGENE"
  try: assocGeneList = genesDict[methGene]
  except Exception: assocGeneList = []
  try: methPos = int(ll[2])
  except Exception: continue
  chrom = "chr"+ll[6]
  try: meanMsc = float(ll[17])
  except Exception: continue
  try: meanMscIfn = float(ll[26])
  except Exception: continue
  meanDiff = meanMscIfn - meanMsc
  try: pvalue = float(ll[30])
  except Exception: continue

  # Filter1
  if(pvalue >= pValueCutoff): continue

  # Filter2
  if(meanDiff > hyperCutoff): # Hyper
    outputFileHyperMeth.write("\t".join([chrom,str(methPos-cpgExtend),str(methPos+cpgExtend)])+"\n")
    for vec in assocGeneList: outputFileHyperMethGenes.write("\t".join([chrom,vec[0],vec[1]])+"\n")
  elif(meanDiff < hypoCutoff): # Hypo
    outputFileHypoMeth.write("\t".join([chrom,str(methPos-cpgExtend),str(methPos+cpgExtend)])+"\n")
    for vec in assocGeneList: outputFileHypoMethGenes.write("\t".join([chrom,vec[0],vec[1]])+"\n")

methFile.close()
outputFileHyperMeth.close()
outputFileHypoMeth.close()
outputFileHyperMethGenes.close()
outputFileHypoMethGenes.close()

# Sorting output files
outputFileNameHyperMethSort = outLoc+"hyper_meth_sort.bed"
outputFileNameHypoMethSort = outLoc+"hypo_meth_sort.bed"
outputFileNameHyperMethGenesSort = outLoc+"hyper_meth_genes_sort.bed"
outputFileNameHypoMethGenesSort = outLoc+"hypo_meth_genes_sort.bed"
toRemove.append(outputFileNameHyperMethSort); toRemove.append(outputFileNameHypoMethSort)
toRemove.append(outputFileNameHyperMethGenesSort); toRemove.append(outputFileNameHypoMethGenesSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameHyperMethTemp+" | uniq > "+outputFileNameHyperMethSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameHypoMethTemp+" | uniq > "+outputFileNameHypoMethSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameHyperMethGenesTemp+" | uniq > "+outputFileNameHyperMethGenesSort)
os.system("sort -k1,1 -k2,2n "+outputFileNameHypoMethGenesTemp+" | uniq > "+outputFileNameHypoMethGenesSort)

# Intersecting with open chromatin region
os.system("intersectBed -a "+outputFileNameHyperMethSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameHyperMeth)
os.system("intersectBed -a "+outputFileNameHypoMethSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameHypoMeth)
os.system("intersectBed -a "+outputFileNameHyperMethGenesSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameHyperMethGenes)
os.system("intersectBed -a "+outputFileNameHypoMethGenesSort+" -b "+dnasePeakFileName+" -wa -u > "+outputFileNameHypoMethGenes)

# Removing files
for e in toRemove: os.system("rm "+e)


