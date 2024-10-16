"""
Read in results from IBDmix
Author: Sara Mathieson
Date: 10/16/24
"""

# python imports
import sys

# our imports
sys.path.append('../gnomad/')
import mask
from region_classes import Individual, Region

################################################################################
# GLOBALS
################################################################################

COMPARE_PATH = "/homes/smathieson/Documents/arg-ml/"
POP = "CEU"
ibdmix_pred_filename = COMPARE_PATH + "IBDmix/Neanderthal_sequence_in_1000genome.50kb_noheader_chr_hg38.txt"
BED_FILE = "/homes/smathieson/Documents/arg-ml/20160622.allChr.mask.bed"

mask_dict = mask.read_mask(BED_FILE)

################################################################################
# HELPERS
################################################################################

def read_ibdmix(target_chrom, pred_filename, target_pop, mask_dict):
    pred_file = open(pred_filename,'r')
    ibdmix_results = {}

    for line in pred_file:
        tokens = line.split()
        chrom = tokens[0][3:]
        pop = tokens[6]
        if chrom == target_chrom and pop == target_pop:
            
            id = tokens[8]
            start = int(tokens[1])
            end = int(tokens[2])
            prob = float(tokens[3]) # LOD not prob...

            region = Region(chrom, start, end, prob)

            if mask.Region(chrom, start, end).inside_mask(mask_dict):

                if id not in ibdmix_results:
                    ibdmix_results[id] = Individual(target_chrom)
            
                ibdmix_results[id].add_region(region)

    pred_file.close()

    # merge haplotypes
    for id in ibdmix_results:
        ibdmix_results[id].collapse_haps() # merge haplotypes

    return ibdmix_results

################################################################################
# MAIN
################################################################################

def ibdmix_one_chrom(CHR):
    ibdmix_results = read_ibdmix(CHR, ibdmix_pred_filename, POP, mask_dict)
    return ibdmix_results

if __name__ == "__main__":
    for chr_int in range(1,23):
        CHR = str(chr_int)
        print(CHR)
        ibdmix_results = ibdmix_one_chrom(CHR)