"""
Read in results from Skov's HMM method
Author: Sara Mathieson
Date: 10/16/24
"""

# python imports
import sys

# our imports
sys.path.append('../gnomad/')
import comparison_helpers
import mask
from region_classes import Individual, Region, read_chrom_lengths

################################################################################
# GLOBALS
################################################################################

COMPARE_PATH = "/homes/smathieson/Documents/arg-ml/"
POP = "CEU"
skov_pred_filename = COMPARE_PATH + "skovHMM/CEU_hg38.txt"
BED_FILE = "/homes/smathieson/Documents/arg-ml/20160622.allChr.mask.bed"

mask_dict = mask.read_mask(BED_FILE)

################################################################################
# HELPERS
################################################################################

def read_skov(target_chrom, pred_filename, target_pop, mask_dict):
    pred_file = open(pred_filename,'r')
    header = pred_file.readline()
    skov_results = {}
    assert target_pop == "CEU" # we only have CEU

    for line in pred_file:
        tokens = line.split()
        chrom = tokens[1][3:]
        if chrom == target_chrom:
            id = tokens[0].split(".")[0]
            start = int(tokens[2])
            end = int(tokens[3])
            prob = float(tokens[4])
            if prob >= 0.8: # their recommended threshold
                region = Region(chrom, start, end, prob)

                if mask.Region(chrom, start, end).inside_mask(mask_dict):

                    if id not in skov_results:
                        skov_results[id] = Individual(target_chrom)
                
                    skov_results[id].add_region(region)
    
    pred_file.close()

    for id in skov_results: # merge haplotypes
        skov_results[id].collapse_haps()

    return skov_results

################################################################################
# MAIN
################################################################################

def skovhmm_one_chrom(CHR):
    skovhmm_results = read_skov(CHR, skov_pred_filename, POP, mask_dict)
    return skovhmm_results

if __name__ == "__main__":
    CHROM_DICT = read_chrom_lengths()

    for chr_int in range(1,23):
        CHR = str(chr_int)
        print("\n---------- starting CHR",CHR,"----------\n")
        results = skovhmm_one_chrom(CHR)

        # optionally compare to nea/san
        comparison_helpers.compare_outgroup_avg(results, "NEA")
        comparison_helpers.compare_outgroup_avg(results, "SAN")

        print("random control")
        random_results = comparison_helpers.make_random_regions(results, CHROM_DICT[CHR])
        comparison_helpers.compare_outgroup_avg(random_results, "NEA")
        comparison_helpers.compare_outgroup_avg(random_results, "SAN")