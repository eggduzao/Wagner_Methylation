
# Import
import os

# Parameters
organism="--organism mm9"
promoter_length="--promoter-length 1000"
maximum_association_length="--maximum-association-length 50000"
multiple_test_alpha="--multiple-test-alpha 0.05"
processes="--processes 1"
print_thresh="--print-thresh 1.0"
bigbed="--bigbed"

###################################################################################################
# PEAKS
###################################################################################################

# Parameters
inputMatrix="/home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/matrix/ME_footprints_peaks.txt"
matchLoc="/home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/results/peaks/"
output_location="--output-location /home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/results/peaks/"

# Execution
clusterCommand = "rgt-motifanalysis --enrichment "
clusterCommand += organism+" "+promoter_length+" "+maximum_association_length+" "+multiple_test_alpha+" "
clusterCommand += processes+" "+output_location+" "+print_thresh+" "+bigbed+" "+inputMatrix+" "+matchLoc
os.system(clusterCommand)

###################################################################################################
# PEAKS + GENES
###################################################################################################

# Parameters
inputMatrix="/home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/matrix/ME_footprints_peaks_genes.txt"
matchLoc="/home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/results/peaks+genes/"
output_location="--output-location /home/egg/Projects/Ash2l/Results/ash2l_mouse/motif_analysis/results/peaks+genes/"

# Execution
clusterCommand = "rgt-motifanalysis --enrichment "
clusterCommand += organism+" "+promoter_length+" "+maximum_association_length+" "+multiple_test_alpha+" "
clusterCommand += processes+" "+output_location+" "+print_thresh+" "+bigbed+" "+inputMatrix+" "+matchLoc
os.system(clusterCommand)


