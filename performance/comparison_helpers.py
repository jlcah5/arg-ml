"""
From a dictionary of regions (for each indvidual), compare 
to the Vindija neanderthal
Author: Sara Mathieson
Date: 10/17/24
"""

# python imports
import numpy as np
import random
from subprocess import Popen, PIPE

# our imports
from region_classes import Individual, Region

################################################################################
# GLOBALS
################################################################################

# order: Altai, Vindija, Denisovan, TODO
NEA_PATH = "/bigdata/smathieson/1000g-share/archaic/"
MOD_PATH = "/homes/smathieson/Documents/arg-ml/phased_haplotypes_v2/"
NEA = ['Vindija33.19']
SAN = ['HGDP00992', 'HGDP01029', 'LP6005441-DNA_A11', 'LP6005441-DNA_B11',
       'LP6005443-DNA_G08', 'SS6004473']

################################################################################
# HELPERS
################################################################################

def parse_bcftools_output(str_output):
    """Return genotype dictionary of all SNPs in the region, for one indv
    format: {<snp pos>: <gt one of 0,1,2>}"""

    lines = str_output.split("\n")
    all_snps = {}
    for l in lines:
        tokens = l.split("\t")
        if len(tokens) > 1:
            gt = tokens[-1]
            if ("|" not in gt) and ("/" not in gt):
                # no samples
                return all_snps
            #print(gt)
            if gt[0] != "." and gt[2] != ".": # not missing
                gt_sum = int(gt[0]) + int(gt[2]) # i.e. "0/1 -> 1"
                snp = str(tokens[1])
                all_snps[snp] = gt_sum
    return all_snps

def bcftools_one_sample(filename, indv_id, CHR, start, end):
    # read vcf file
    cmd = "bcftools view --no-header --force-samples -s " + indv_id + " -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + filename
    process = Popen(cmd, shell=True, stdout=PIPE)
    output, err = process.communicate()
    output = output.decode("utf-8")
    snps = parse_bcftools_output(output)
    return snps

def get_identity(snpsA, snpsB):

    overlap = set(snpsA.keys()).intersection(set(snpsB.keys()))
    #print(len(overlap))

    # identity (not percent!)
    identity = 0
    num_skip = 0
    for snp in overlap:
        nea_gt = snpsA[snp]
        ceu_gt = snpsB[snp]
        #print(nea_gt, ceu_gt)

        # 0&0 add 1, 1&1 add 1, 2&2 add 1
        # 0&1 add 0.5, 1&2 add 0.5
        # 0&2 add 0
        diff = abs(nea_gt - ceu_gt)
        if diff == 0:
            identity += 1
        elif diff == 1:
            identity += 0.5
        elif diff == 2:
            identity += 0
        else:
            #print("unknown diff", nea_gt, ceu_gt)
            num_skip += 1
    
    #identity = identity/(len(overlap)-num_skip)
    #print("percent identity", identity)
    return identity, len(overlap)-num_skip
    #input('enter')

def compare_other(indv_id, CHR, start, end, flag):
    # one individual, one region
    nea_filename = NEA_PATH + "individuals_highcov." + CHR + ".bcf"
    hgdp_filename = MOD_PATH + "hgdp1kgp_chr" + CHR + ".filtered.SNV_INDEL.phased.shapeit5.bcf"

    # read nea
    if flag == "NEA":
        other_snps = bcftools_one_sample(nea_filename, NEA[0], CHR, start, end)

    # read san
    elif flag == "SAN":
        other_snps = bcftools_one_sample(hgdp_filename, SAN[0], CHR, start, end)

    else:
        print("unknown flag", flag)
        sys.exit()

    # read ceu
    ceu_snps = bcftools_one_sample(hgdp_filename, indv_id, CHR, start, end)

    return get_identity(other_snps, ceu_snps)

def compare_outgroup_avg(region_dict, flag):

    avg_all = 0
    num_skip = 0
    for id in region_dict:
        indv_avg = 0
        indv_snp = 0
        for region in region_dict[id].region_lst:
            identity, num_snps = compare_other(id, region.chrom, region.start, region.end, flag)
            indv_avg += identity
            indv_snp += num_snps
        
        if indv_snp == 0: # indv not in file
            num_skip += 1
        else:
            indv_avg = indv_avg/indv_snp
            #print(id, indv_avg)
            avg_all += indv_avg
        #input('enter')
    
    avg_all = avg_all/(len(region_dict)-num_skip)
    print(flag, "overall avg percent identity", avg_all)

def make_random_regions(region_dict, chrom_len):
    """From a given region_dict, create new random regions of the same lengths
    """
    random_regions = {}
    for id in region_dict:

        rand_region_lst = []
        for region in region_dict[id].region_lst:
            new_start = random.randrange(1,chrom_len-region.end)
            new_end = new_start + (region.end - region.start)
            new_region = Region(region.chrom, new_start, new_end, region.prob)
            rand_region_lst.append(new_region)
        
        new_indv = Individual(region.chrom)
        new_indv.region_lst = rand_region_lst
        random_regions[id] = new_indv
    
    return random_regions


#def compare_nea(CHR, region_dict):

################################################################################
# MAIN
################################################################################

#if __name__ == "__main__":
    # testing
    # length of chrom
    #chrom_dict = read_chrom_lengths()
    #read_nea("1", chrom_dict)
