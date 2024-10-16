"""
Read in results from our method (EGRM-based)
Author: Sara Mathieson
Date: 10/16/24
"""

# our imports
from region_classes import Individual, Region

################################################################################
# GLOBALS
################################################################################

PRED_PATH = "/homes/smathieson/Documents/arg-ml/output/"
ID_PATH = "/homes/smathieson/GIT/arg-ml/gnomad/gnomad_subpops/"
POP = "CEU"

# our pred file and ID file
egrm_pred_path = PRED_PATH + POP + "/"
egrm_id_filename = ID_PATH + POP.lower() + "_yri.txt"

THRESH = 0.997

################################################################################
# HELPERS
################################################################################

def read_egrm(target_chrom, pred_filename, id_filename): # don't need mask
    # first is target, second is outgroup
    id_lst = open(id_filename, 'r').read().split()
    id_lst = id_lst[:len(id_lst)//2] # just take target
    
    pred_file = open(pred_filename, 'r')
    egrm_results = {}
    for line in pred_file:
        tokens = line.split()
        prob = float(tokens[3])

        # only take high prob regions for us
        if prob >= THRESH:
            hap_id = int(tokens[4])
            id = id_lst[hap_id//2]
            chrom = tokens[0]
            assert target_chrom == chrom
            start = int(tokens[1])
            end = int(tokens[2])

            region = Region(chrom, start, end, prob)

            if id not in egrm_results:
                egrm_results[id] = Individual(target_chrom)
            
            egrm_results[id].add_region(region)

    pred_file.close()

    # merge consecutive regions and haplotypes
    for id in egrm_results:
        indv = egrm_results[id]
        indv.merge_consecutive()
        indv.collapse_haps() # merge haplotypes
    return egrm_results

################################################################################
# MAIN
################################################################################

def egrm_one_chrom(CHR):
    egrm_pred_filename = egrm_pred_path + POP + "_chr" + CHR + ".pred"
    egrm_results = read_egrm(CHR, egrm_pred_filename, egrm_id_filename)
    return egrm_results

if __name__ == "__main__":

    for chr_int in range(1,23):
        CHR = str(chr_int)
        egrm_results = egrm_one_chrom(CHR)
        