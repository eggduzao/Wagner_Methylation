
# Import
import os

# Parameters
organism="--organism hg19"
fpr="--fpr 0.0001"
precision="--precision 10000"
pseudocounts="--pseudocounts 0.1"
rand_proportion="--rand-proportion 10"


# Execution - Genes
loc="/Users/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/"
inputList = [loc+"genes_down.bed", loc+"genes_up.bed", loc+"hyper_meth_genes.bed", loc+"hypo_meth_genes.bed"]
output_location="--output-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/genes/match/"
command = "rgt-motifanalysis --matching "+" ".join([organism,fpr,precision,pseudocounts,rand_proportion,output_location]+inputList)
#os.system(command)

# Executuon - Methyl
loc="/Users/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/"
inputList = [loc+"hyper_meth.bed", loc+"hypo_meth.bed"]
output_location="--output-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/methylation/match/"
command = "rgt-motifanalysis --matching "+" ".join([organism,fpr,precision,pseudocounts,rand_proportion,output_location]+inputList)
os.system(command)


