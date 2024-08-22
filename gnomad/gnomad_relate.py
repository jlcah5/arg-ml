"""
Description: read combined hgdp+1kgp data, split into populations, infer the
    ARG, and calculate eGRM on sliding windows (window size set to 50kb, step
    size to 10kb)
Usage: python3 gnomad_relate.py IN_FOLDER BED_FILE RECO_FOLDER OUT_FOLDER
Author: Jordan Cahoon, Sara Mathieson
Date: 7/17/24
"""

# python imports
import numpy as np
from subprocess import Popen, PIPE
import sys

# our imports
import mask

################################################################################
# GLOBALS
################################################################################

IN_FOLDER = sys.argv[1]
POP = IN_FOLDER[-3:]
BED_FILE = sys.argv[2]
RECO_FOLDER = sys.argv[3] + "/" + POP
OUT_FOLDER = sys.argv[4] + "/" + POP
START_CHR = 1
END_CHR = 22 # inclusive

# start or end part way through
if len(sys.argv) == 6:
    START_CHR = int(sys.argv[5])
elif len(sys.argv) == 7:
    START_CHR = int(sys.argv[5])
    END_CHR = int(sys.argv[6])

# set up directory
mkdir = Popen("mkdir " + OUT_FOLDER, shell=True, stdout=PIPE)
mkdir.communicate()

WINDOW = 50000 # 50kb
STEP   = 10000 # 10kb
MIN_SNPS = 50 # min SNPs per 50kb region

################################################################################
# HELPERS
################################################################################

def read_chrom_lengths():
    arr = np.loadtxt("hg38_chrom_lengths.tsv", dtype='int', delimiter="\t", skiprows=1)
    chrom_dict = {}
    for chr in range(1,23):
        assert arr[chr-1][0] == chr
        chrom_dict[str(chr)] = arr[chr-1][1]
    return chrom_dict

################################################################################
# MAIN
################################################################################

def main():

    # length of chrom
    chrom_dict = read_chrom_lengths()

    # accessibility mask
    mask_dict = mask.read_mask(BED_FILE)
    kept = 0
    total = 0

    for chr_int in range(START_CHR,END_CHR+1):
        CHR = str(chr_int)
        chrom_length = chrom_dict[CHR]

        # go through each region
        start = 0
        end = WINDOW
        
        while end <= chrom_length:
            total += 1

            cmd = "bcftools view --no-header -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + IN_FOLDER + "/" + POP + "_chr" + CHR + ".vcf.gz | wc -l"
            process = Popen(cmd, shell=True, stdout=PIPE)
            output, err = process.communicate()
            num_snps = int(output.decode("utf-8"))

            # create region to determine accessibility
            region = mask.Region(CHR, start, end)
            
            # if we have enough SNPs and inside accessibility mask
            if num_snps >= MIN_SNPS and region.inside_mask(mask_dict):
                kept += 1

                print("num SNPs", num_snps)
                prefix = POP + "_chr" + CHR + "_" + str(start) + "_" + str(end)

                # extract SNPs
                bcftools_cmd = "bcftools view -r chr" + CHR + ":" + str(start) + "-" + str(end) + " -Oz -o " + prefix + ".vcf.gz " + IN_FOLDER + "/" + POP + "_chr" + CHR + ".vcf.gz"
                process = Popen(bcftools_cmd, shell=True, stdout=PIPE)
                process.communicate() # wait to finish

                # convert from vcf.gz to haps and sample
                relate = "RelateFileFormats --mode ConvertFromVcf --haps " + prefix + ".haps --sample " + prefix + ".sample -i " + prefix
                process = Popen(relate, shell=True, stdout=PIPE)
                process.communicate() # wait to finish

                # infer tree (biallelic snps retained with bcftools)
                # use population specific genetic map from Spence and Song
                reco_file = RECO_FOLDER + "/" + POP + "_recombination_map_hapmap_format_hg38_chr_" + CHR + ".txt"
                relate = "Relate --mode All -m 1.25e-8 -N 20000 --haps " + prefix + ".haps --sample " + prefix + ".sample --map " + reco_file + " --seed 1 -o " + prefix
                process = Popen(relate, shell=True, stdout=PIPE)
                process.communicate() # wait to finish
                
                # convert to tree sequence
                relate = "RelateFileFormats --mode ConvertToTreeSequence -i " + prefix + " -o " + prefix + ".infer"
                process = Popen(relate, shell=True, stdout=PIPE)
                process.communicate() # wait to finish

                # calculate GRM on entire region
                egrm = "trees2egrm --output-format numpy " + prefix + ".infer.trees --c --haploid --output " + OUT_FOLDER + "/" + prefix
                process = Popen(egrm, shell=True, stdout=PIPE)
                process.communicate() # wait to finish
                
                # clean 
                clean = "rm -rf " + prefix + "*"
                process = Popen(clean, shell=True, stdout=PIPE)
                process.communicate() # wait to finish

            start += STEP
            end += STEP

    print("frac kept", kept/total)

if __name__ == "__main__":
    main()

