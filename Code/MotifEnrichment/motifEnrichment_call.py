
# Import
import os

# Parameters
organism="--organism hg19"
promoter_length="--promoter-length 1000"
maximum_association_length="--maximum-association-length 50000"
multiple_test_alpha="--multiple-test-alpha 0.05"
print_thresh="--print-thresh 1.0"
no_copy_logos="--no-copy-logos"

# Execution Genes
backgroundBedFile = "/Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/genes/match/random_regions.bed"
loc="/Users/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/"
inputList = [loc+"genes_down.bed", loc+"genes_up.bed", loc+"hyper_meth_genes.bed", loc+"hypo_meth_genes.bed"]
matching_location="--matching-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/genes/match/"
output_location="--output-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/genes/enrichment/"
command = "rgt-motifanalysis --enrichment "+" ".join([organism,promoter_length,maximum_association_length,multiple_test_alpha,matching_location,output_location,print_thresh,no_copy_logos,backgroundBedFile]+inputList)
os.system(command)

# Execution Methyl
backgroundBedFile = "/Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/methylation/match/random_regions.bed"
loc="/Users/egg/Projects/Wagner_Methylation/Results/Parsed_Regions/"
inputList = [loc+"hyper_meth.bed", loc+"hypo_meth.bed"]
matching_location="--matching-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/methylation/match/"
output_location="--output-location /Users/egg/Projects/Wagner_Methylation/Results/MotifEnrichment/methylation/enrichment/"
command = "rgt-motifanalysis --enrichment "+" ".join([organism,promoter_length,maximum_association_length,multiple_test_alpha,matching_location,output_location,print_thresh,no_copy_logos,backgroundBedFile]+inputList)
os.system(command)


