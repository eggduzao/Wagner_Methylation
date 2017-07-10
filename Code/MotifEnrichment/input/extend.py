import os
os.system("slopBed -i result.bed -g /hpcwork/izkf/projects/TfbsPrediction/Data/HG19/hg19.chrom.sizes.filtered -b 25 > temp.bed")
os.system("cut -f 1,2,3,4 temp.bed > result_extended.bed")
os.system("rm temp.bed")
