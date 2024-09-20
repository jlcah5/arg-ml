"""
Compare num SNPs between real and simulated data.
Based on simulator-Distribution-sw.py and simulator-Distribution-sw.sh
Author: Sara Mathieson
Date: 9/18/24
"""

import demes
#import msprime
import numpy as np
import random
from subprocess import Popen, PIPE
import sys

# our imports
import mask

# intializing variables for sims
length=200000 # like a chromosome
num_models = 200 # number of demographies from unified demographic model
num_sims = 2
model_number = random.randrange(0,num_models) # choose a random one
print("model number", model_number)
model_type = "growth"
window_size=50000 # 50kb
step=10000 # 10kb
n=56 # number of individuals

# initializing variables for real
WINDOW = 50000 # 50kb
STEP   = 10000 # 10kb
MIN_SNPS = 50 # min SNPs per 50kb region
IN_FOLDER = sys.argv[1]
BED_FILE = sys.argv[2]
POP = IN_FOLDER[-3:] # CEU
print(POP)

def read_chrom_lengths():
    arr = np.loadtxt("hg38_chrom_lengths.tsv", dtype='int', delimiter="\t", skiprows=1)
    chrom_dict = {}
    for chr in range(1,23):
        assert arr[chr-1][0] == chr
        chrom_dict[str(chr)] = arr[chr-1][1]
    return chrom_dict

def test_sims():
    # loading demographic model
    graph = demes.load(f'../simulation/demographic_models/distribution_{model_type}/yaml/model_{model_number}.yaml')

    demography = msprime.Demography.from_demes(graph)
    demography.sort_events()

    for i in range(num_sims): 
        # kept outgroup at n to compare with real data
        ts = msprime.sim_ancestry(samples={'ANC':0, "SRC": 0 , "TAR": n, "OUT": n}, demography=demography, sequence_length=length, ploidy=2, recombination_rate = 1.25e-8, record_migrations=True)
        ts = msprime.sim_mutations(ts, rate = 1.29e-8, discrete_genome = False)

        for start in range(0,length-window_size+step,step): # create files for each 50kb window
            ts_new = ts.keep_intervals(np.array([[start, start+window_size]]), simplify=False).trim()
            print(ts_new.num_mutations)
            
def test_real():

    # length of chrom
    chrom_dict = read_chrom_lengths()

    # accessibility mask
    mask_dict = mask.read_mask(BED_FILE)

    for chr_int in range(1,23):
        CHR = str(chr_int)
        chrom_length = chrom_dict[CHR]

        # go through each region
        start = random.randrange(0, chrom_length)
        end = start + WINDOW
        print(CHR, start, end)
        
        #while end <= chrom_length:
            #total += 1

        cmd = "bcftools view --no-header -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + IN_FOLDER + "/" + POP + "_chr" + CHR + ".vcf.gz | wc -l"
        process = Popen(cmd, shell=True, stdout=PIPE)
        output, err = process.communicate()
        num_snps = int(output.decode("utf-8"))

        # create region to determine accessibility
        region = mask.Region(CHR, start, end)
        
        # if we have enough SNPs and inside accessibility mask
        if num_snps >= MIN_SNPS and region.inside_mask(mask_dict):
            print(num_snps)
            input('enter')

if __name__ == "__main__":
    #test_sims()
    test_real()

