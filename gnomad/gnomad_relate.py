"""
Description: read combined hgdp+1kgp data, split into populations, infer the
    ARG, and calculate eGRM on sliding windows (window size set to 50kb, step
    size to 10kb)
Usage: python3 gnomad_relate.py IN_FOLDER OUT_FOLDER
Author: Jordan Cahoon, Sara Mathieson
Date: 6/12/24

TODO:
- adjust relate simulation
- create tarball for input into model
"""

# python imports
from subprocess import Popen, PIPE
import sys

# our imports
sys.path.insert(1, "../../pg-gan")
import real_data_random

################################################################################
# GLOBALS
################################################################################

# TODO expand to all chroms and pops
CHR = "3"
POP = "ACB"

IN_FOLDER = sys.argv[1]
BED_FILE = sys.argv[2]
OUT_FOLDER = sys.argv[3] + "/" + POP

WINDOW = 50000 # 50kb
STEP   = 10000 # 10kb
MIN_SNPS = 50 # min SNPs per 50kb region

################################################################################
# HELPERS
################################################################################

def read_chrom_lengths():
    arr = np.loadtxt("hg38_chrom_lengths.txt", dtype='int', delimiter="\t", skip_rows=1)
    chrom_dict = {}
    for chr in range(1,23):
        assert arr[chr-1][0] = chr
        chrom_dict[str(chr)] = arr[chr-1][1]
    print(chrom_dict)
    return chrom_dict

################################################################################
# MAIN
################################################################################

def main():

    # length of chrom
    chrom_dict = read_chrom_lengths():
    chrom_length = chrom_dict[CHR]

    # accessibility mask
    mask_dict = real_data_random.read_mask(BED_FILE)

    # go through each region
    start = 0
    end = WINDOW
    kept = 0
    while end <= chrom_length:

        pop_file = "gnomad_subpops/" + POP.lower() + ".txt"
        cmd = "bcftools view --no-header -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + IN_FOLDER + "/" + POP + "/" + POP + "_chr" + CHR + ".vcf.gz | wc -l"
        process = Popen(cmd, shell=True, stdout=PIPE)
        output, err = process.communicate()
        num_snps = int(output.decode("utf-8"))

        # create region to determine accessibility
        region = real_data_random.Region(CHR, start, end)
        
        # if we have enough SNPs and inside accessibility mask
        if num_snps >= MIN_SNPS and region.inside_mask(mask_dict):
            print("num SNPs", num_snps)
            kept += 1
            prefix = POP + "_chr" + CHR + "_" + str(start) + "_" + str(end)

            # extract SNPs
            bcftools_cmd = "bcftools view -r chr" + CHR + ":" + str(start) + "-" + str(end) + " -Oz -o " + prefix + ".vcf.gz " + IN_FOLDER + "/" + POP + "/" + POP + "_chr" + CHR + ".vcf.gz"
            process = Popen(bcftools_cmd, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # convert from vcf.gz to haps and sample
            relate = "RelateFileFormats --mode ConvertFromVcf --haps " + prefix + ".haps --sample " + prefix + ".sample -i " + prefix
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # infer tree (biallelic snps retained with bcftools))
            relate = "Relate --mode All -m 1.25e-8 -N 20000 --haps " + prefix + ".haps --sample " + prefix + ".sample --map ../simulation/genetic_map.txt --seed 1 -o " + prefix
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish
            
            # convert to tree sequence
            relate = "RelateFileFormats --mode ConvertToTreeSequence -i " + prefix + " -o " + prefix + ".infer"
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # calculate GRM on entire region
            egrm = "trees2egrm --output-format numpy " + prefix + ".infer.trees --c --haploid --output " + OUT_FOLDER + "/" + prefix #" --left ${win_start} --right ${win_end}
            process = Popen(egrm, shell=True, stdout=PIPE)
            process.communicate() # wait to finish
            
            # clean 
            clean = "rm -rf " + prefix + "*"
            print(clean)
            process = Popen(clean, shell=True, stdout=PIPE)
            process.communicate() # wait to finish
            input('enter')

        start += STEP
        end += STEP

    print("frac kept", kept/1000)

# TODO: tarball of PATH_TO_OUTPUT/*.npy files

if __name__ == "__main__":
    main()

