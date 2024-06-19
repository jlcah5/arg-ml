"""
Description: separate VCFs by population and keep only biallelic SNPs
Usage: python3 prep_vcf.py IN_FOLDER
Author: Sara Mathieson
Date: 6/14/24
"""

# python imports
import os
from subprocess import Popen, PIPE
import sys

################################################################################
# GLOBALS
################################################################################

IN_FOLDER = sys.argv[1]

################################################################################
# MAIN
################################################################################

def main():

    pop_file_lst = os.listdir("gnomad_subpops")
    pop_file_lst.sort()

    # go through each population
    for pop_file in pop_file_lst:
        if "_yri.txt" in pop_file:
            pop = pop_file[:3].upper()
            pop_file = "gnomad_subpops/" + pop_file

            # folder for each pop
            out_folder = IN_FOLDER + "/" + pop
            mkdir = "mkdir " + out_folder
            process = Popen(mkdir, shell=True, stdout=PIPE)
            process.communicate()

            # go through each chromosome
            for i in range(1,23):
                chr = str(i)
                in_file = IN_FOLDER + "/hgdp1kgp_chr" + chr + ".filtered.SNV_INDEL.phased.shapeit5.bcf"
                out_file = out_folder + "/" + pop + "_chr" + chr + ".vcf.gz"
                
                view_norm = "bcftools view -S " + pop_file + " --force-samples --min-ac 1:minor -m2 -M2 -v snps -Ou " + in_file + " | bcftools norm -d snps -Oz -o " + out_file
                print(view_norm)
                process = Popen(view_norm, shell=True, stdout=PIPE)
                process.communicate()

                # index
                index = "bcftools index " + out_folder + "/" + pop + "_chr" + chr + ".vcf.gz"
                print(index)
                process = Popen(index, shell=True, stdout=PIPE)
                process.communicate()

if __name__ == "__main__":
    main()