"""
Compare num SNPs between real and simulated data.
Based on simulator-Distribution-sw.py and simulator-Distribution-sw.sh
Author: Sara Mathieson
Date: 9/18/24
"""

import demes
import matplotlib.pyplot as plt
import msprime
import numpy as np
import random
import seaborn as sns
from subprocess import Popen, PIPE
import sys

# our imports
import mask

num_regions = 5000
WINDOW = 50000 # 50kb
STEP   = 10000 # 10kb

# intializing variables for sims
length=200000 # like a chromosome
num_models = 200 # number of demographies from unified demographic model
num_sims = int(num_regions/(length/STEP))
model_type = "growth"
N=56 # number of individuals

# initializing variables for real
MIN_SNPS = 10 # min SNPs per 50kb region TODO change to 50!
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
    
    snp_counts = []
    for i in range(num_sims*2): # add more since some models are not working 
        if i % num_sims*2 == 10:
            print(i, "/", num_sims*2)

        # set up demographic model
        model_number = random.randrange(0,num_models)
        try:
            graph = demes.load(f'../simulation/demographic_models/distribution_{model_type}/yaml/model_{model_number}.yaml')
            demography = msprime.Demography.from_demes(graph)
            demography.sort_events()

            # kept outgroup at n to compare with real data
            ts = msprime.sim_ancestry(samples={'ANC':0, "SRC": 0 , "TAR": N, "OUT": N}, demography=demography, sequence_length=length, ploidy=2, recombination_rate = 1.25e-8, record_migrations=True)
            ts = msprime.sim_mutations(ts, rate = 1.29e-8, discrete_genome = False)

            for start in range(0,length-WINDOW+STEP,STEP): # create files for each 50kb window
                ts_new = ts.keep_intervals(np.array([[start, start+WINDOW]]), simplify=False).trim()
                #print(ts_new.num_mutations)
                snp_counts.append(ts_new.num_mutations)
                #input('enter')
        except:
            #print("model number", model_number, "not working")
            pass

    #np.save('sim_snp_counts.npy', snp_counts)
    return snp_counts
            
def test_real():

    # length of chrom
    chrom_dict = read_chrom_lengths()

    # accessibility mask
    mask_dict = mask.read_mask(BED_FILE)

    snp_counts = []
    for x in range(num_regions):
        if x % 100 == 0:
            print(x , "/", num_regions)
        chr_int = random.randrange(1,23)
        CHR = str(chr_int)
        chrom_length = chrom_dict[CHR]

        # go through each region
        start = random.randrange(0, chrom_length)
        end = start + WINDOW

        cmd = "bcftools view --no-header -r chr" + CHR + ":" + str(start) + "-" + str(end) + " " + IN_FOLDER + "/" + POP + "_chr" + CHR + ".vcf.gz | wc -l"
        process = Popen(cmd, shell=True, stdout=PIPE)
        output, err = process.communicate()
        num_snps = int(output.decode("utf-8"))

        # create region to determine accessibility
        region = mask.Region(CHR, start, end)
        
        # if we have enough SNPs and inside accessibility mask
        if num_snps >= MIN_SNPS and region.inside_mask(mask_dict):
            snp_counts.append(num_snps)

    #np.save('real_snp_counts.npy', snp_counts)
    return snp_counts

def plot_sim_real(sim_snps, real_snps):
    #sim_snps = np.load('sim_snp_counts.npy', allow_pickle=True)
    print(len(sim_snps), sim_snps[:10])
    #real_snps = np.load('real_snp_counts.npy', allow_pickle=True)
    print(len(real_snps), real_snps[:10])
    sns.distplot(sim_snps, label="sim")
    sns.distplot(real_snps, label="real")
    plt.legend()
    plt.savefig("test_num_snps.pdf")

if __name__ == "__main__":
    sim_snps = test_sims()
    real_snps = test_real()

    plot_sim_real(sim_snps, real_snps)

