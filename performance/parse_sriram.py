"""
Read in results from Sriram's CRF method
Author: Sara Mathieson
Date: 10/16/24
"""

# python imports
import gzip
from subprocess import Popen, PIPE
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
BED_FILE = "/homes/smathieson/Documents/arg-ml/20160622.allChr.mask.bed"
CHAIN_FILE = "/homes/smathieson/Programs/hg19ToHg38.over.chain.gz"

# accessibility mask
mask_dict = mask.read_mask(BED_FILE)

# files for sriram method
sriram_pred_path = COMPARE_PATH + "sriram/summaries.release/" + POP + ".hapmap/summaries/haplotypes/"
sriram_id_filename = COMPARE_PATH + "sriram/summaries.release/ids/" + POP + ".ids"

################################################################################
# HELPERS
################################################################################

def liftOver(region_lst):
    # convert from hg19 to hg38 using cumbersome procedure..
    
    # write regions to hg19 file
    hg19_filename = "temp_" + POP + "_hg19.txt"
    hg19_file = open(hg19_filename, 'w')
    for region in region_lst:
        hg19_file.write(" ".join(["chr" + region.chrom, str(region.start), str(region.end), str(region.prob)]) + "\n")
    hg19_file.close()

    # liftOver command
    hg38_filename = "temp_" + POP + "_hg38.txt"
    liftOver_cmd = "liftOver " + hg19_filename + " " + CHAIN_FILE + " " + hg38_filename + " unMapped" + POP
    process = Popen(liftOver_cmd, shell=True, stdout=PIPE)
    process.communicate() # wait to finish

    # read from hg38 file
    hg38_file = open(hg38_filename, 'r')
    new_regions = []
    for line in hg38_file:
        new_chr, new_start, new_end, new_prob = line.split()
        region = Region(new_chr[3:], int(new_start), int(new_end), float(new_prob))
        new_regions.append(region)
    hg38_file.close()
    
    # delete temp_hg19.txt temp_hg38.txt and unMapped
    clean_up = "rm " + hg19_filename + " " + hg38_filename + " unMapped" + POP
    process = Popen(clean_up, shell=True, stdout=PIPE)
    process.communicate() # wait to finish

    return new_regions

def read_sriram(target_chrom, pred_filename, id_filename, mask_dict, CHAIN_FILE):
    # first get 1000g IDs
    id_file = open(id_filename, 'r')
    id_lst = []
    for line in id_file:
        tokens = line.split()
        id_lst.append(tokens[0])
    id_file.close()
    
    pred_file = gzip.open(pred_filename, 'rt') # gz file
    sriram_results = {}
    for line in pred_file:
        if not line.startswith("##"):
            tokens = line.split()
            hap_id = int(tokens[1])
            id = id_lst[hap_id]
            chrom = tokens[0]
            assert target_chrom == chrom
            start = int(tokens[2])
            end = int(tokens[3])
            prob = float(tokens[6])

            region = Region(chrom, start, end, prob)

            if mask.Region(chrom, start, end).inside_mask(mask_dict):

                if id not in sriram_results:
                    sriram_results[id] = Individual(target_chrom)
            
                sriram_results[id].add_region(region)
    
    pred_file.close()

    # convert to hg38 (liftOver) + merge haplotypes
    for id in sriram_results:
        indv = sriram_results[id]
        indv.region_lst = liftOver(indv.region_lst)
        indv.collapse_haps() # merge haplotypes

    return sriram_results

################################################################################
# MAIN
################################################################################

def sriram_one_chrom(CHR):
    sriram_pred_filename = sriram_pred_path + "chr-" + CHR + ".thresh-90.length-0.00.haplotypes.gz"
    sriram_results = read_sriram(CHR, sriram_pred_filename, sriram_id_filename, mask_dict, CHAIN_FILE)
    return sriram_results

if __name__ == "__main__":
    
    for chr_int in range(1,23):
        CHR = str(chr_int)
        sriram_results = sriram_one_chrom(CHR)


        