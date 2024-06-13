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

CHR = "21"
POP = "CEU"

IN_FOLDER = sys.argv[1]
OUT_FOLDER = sys.argv[2] + "/" + POP

WINDOW = 50000 # 50kb
STEP   = 10000 # 10kb

################################################################################
# MAIN
################################################################################

def main():

    # go through each region
    start = 0
    end = WINDOW
    for i in range(1000): # TODO how to determine end

        pop_file = "gnomad_subpops/" + POP.lower() + ".txt"
        cmd = "bcftools view --no-header -S " + pop_file + " --force-samples --min-ac 1:minor -m2 -M2 -v snps -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + IN_FOLDER + "/hgdp1kgp_chr" + CHR + ".filtered.SNV_INDEL.phased.shapeit5.bcf | wc -l"
        process = Popen(cmd, shell=True, stdout=PIPE)
        output, err = process.communicate()
        num_snps = int(output.decode("utf-8"))
        
        # no SNPs
        if num_snps < 20:
            pass
        else:
            print("num SNPs", num_snps)
            prefix = POP + "_chr" + CHR + "_" + str(start) + "_" + str(end)

            # extract SNPs
            bcftools_cmd = "bcftools view -S " + pop_file + " --force-samples --min-ac 1:minor -m2 -M2 -v snps -r chr" + CHR + ":" + str(start) + "-" + str(end) + " -Oz -o " + prefix + ".vcf.gz " + IN_FOLDER + "/hgdp1kgp_chr" + CHR + ".filtered.SNV_INDEL.phased.shapeit5.bcf"
            process = Popen(bcftools_cmd, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # convert from vcf.gz to haps and sample
            relate = "RelateFileFormats --mode ConvertFromVcf --haps " + prefix + ".haps --sample " + prefix + ".sample -i " + prefix
            #print(relate)
            #input('enter')
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # remove non biallelic snps (not needed, done with bcftools)
            #relate = "RelateFileFormats --mode RemoveNonBiallelicSNPs --haps " + prefix + ".haps -o " + prefix + ".clean"
            #Popen(relate, shell=True, stdout=PIPE)

            # infer tree
            relate = "Relate --mode All -m 1.25e-8 -N 20000 --haps " + prefix + ".haps --sample " + prefix + ".sample --map ../simulation/genetic_map.txt --seed 1 -o " + prefix
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish
            
            # convert to tree sequence
            relate = "RelateFileFormats --mode ConvertToTreeSequence -i " + prefix + " -o " + prefix + ".infer"
            process = Popen(relate, shell=True, stdout=PIPE)
            process.communicate() # wait to finish

            # calculate GRM on sliding window
            #for win_start in $( eval echo {${start}..${end}..10000} )
            #do
            #let win_end=$win_start+50000
            #trees2egrm --output-format numpy ${prefix}.infer.trees --c --haploid --output PATH_TO_OUTPUT/$prefix_${win_start}_${win_end} --left ${win_start} --right ${win_end}
            
            #done
            # clean 
            #rm -rf ${prefix}.*'''
            input('enter')

        start += STEP
        end += STEP

# to-do: tarball of PATH_TO_OUTPUT/*.npy files

if __name__ == "__main__":
    main()

