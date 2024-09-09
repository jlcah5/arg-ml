
# python imports
import gzip
import numpy as np
import random
from scipy.stats import norm
from subprocess import Popen, PIPE
import sys

################################################################################
# GLOBALS
################################################################################

# command line args
CHAIN_FILE = sys.argv[1] #"/homes/smathieson/Programs/hg19ToHg38.over.chain.gz"
PRED_PATH = sys.argv[2] #"/homes/smathieson/Documents/arg-ml/output/"
ID_PATH = sys.argv[3] #"/homes/smathieson/GIT/arg-ml/gnomad/gnomad_subpops/"
COMPARE_PATH = sys.argv[4] #/homes/smathieson/Documents/arg-ml/
POP = sys.argv[5]

def read_chrom_lengths():
    arr = np.loadtxt("../gnomad/hg38_chrom_lengths.tsv", dtype='int', delimiter="\t", skiprows=1)
    chrom_dict = {}
    for chr in range(1,23):
        assert arr[chr-1][0] == chr
        chrom_dict[str(chr)] = arr[chr-1][1]
    return chrom_dict

CHROM_DICT = read_chrom_lengths()

THRESH = 0.998
RAND_TRIALS = 1000
SIG_THRESH = 0.05 # significance threshold (p-value)

# our pred file and ID file
egrm_pred_path = PRED_PATH + POP + "/"
egrm_id_filename = ID_PATH + POP.lower() + "_yri.txt"

# files for other methods
sriram_pred_path = COMPARE_PATH + "sriram/summaries.release/" + POP + ".hapmap/summaries/haplotypes/"
sriram_id_filename = COMPARE_PATH + "sriram/summaries.release/ids/" + POP + ".ids"
ibdmix_pred_filename = COMPARE_PATH + "IBDmix/Neanderthal_sequence_in_1000genome.50kb_noheader_chr_hg38.txt"
skov_pred_path = COMPARE_PATH + "skovHMM/CEU_hg38.txt"

################################################################################
# CLASSES
################################################################################

class Individual:

    def __init__(self, chrom):
        # TODO add population, etc? id number?
        self.chrom = chrom
        self.region_lst = []

    def add_region(self, region):
        if region not in self.region_lst:
            self.region_lst.append(region)
        else:
            # TODO could average the probs for each hap or something
            # rare case so skip for now
            pass
        # TODO are these in chromosomal order?

    def merge_consecutive(self):
        # only needed for our method
        merged_idx_lst = []
        curr_idx = [0]
        for i in range(len(self.region_lst)-1):
            end = self.region_lst[i].end
            start = self.region_lst[i+1].start
            if end == start: # end of region is same as start of next region
                #print(self.region_lst[i], self.region_lst[i+1])
                curr_idx.append(i+1)
            else:
                merged_idx_lst.append(curr_idx)
                curr_idx = [i+1]
        merged_idx_lst.append(curr_idx)

        #print(len(merged_idx_lst), len(self.region_lst))
        assert merged_idx_lst[-1][-1] + 1 == len(self.region_lst)

        # use index groups to redo regions, averaging preds
        merged_regions = []
        for idx_group in merged_idx_lst:
            if len(idx_group) == 1:
                merged_regions.append(self.region_lst[idx_group[0]])
            else:
                first_idx = idx_group[0]
                last_idx = idx_group[-1]
                chr = self.region_lst[first_idx].chrom
                start = self.region_lst[first_idx].start
                end = self.region_lst[last_idx].end
                new_prob = np.mean([region.prob for region in self.region_lst[first_idx:last_idx+1]])
                new_region = Region(chr, start, end, new_prob)
                merged_regions.append(new_region)

        self.region_lst = merged_regions

    def frac_nea(self):
        total = CHROM_DICT[self.chrom]
        base_sum = sum([(x.end - x.start) for x in self.region_lst])
        return base_sum/total
            
    def __str__(self):
        return "indv nea frac: " + str(self.frac_nea())

class Region:

    def __init__(self, chrom, start, end, prob):
        self.chrom = chrom # str
        self.start = start
        self.end = end
        self.prob = prob # prob introgression

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            # don't compare prob
            return self.chrom == other.chrom and self.start == other.start and self.end == other.end
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return " ".join([self.chrom, str(self.start), str(self.end), str(self.prob)])

################################################################################
# HELPERS
################################################################################

def shuffle_indvs(result_dict):
    # shuffle called regions between the individuals as a way to get a meta p-value
    # right now this is all within one chrom
    result_items = list(result_dict.items())

    # get keys and values before and after
    keys, values = zip(*result_items)
    random.shuffle(result_items)
    shuffled_keys, shuffled_values = zip(*result_items)

    # make new dictionary with shuffled keys and same values
    random_dict = dict(zip(shuffled_keys, values))
    return random_dict

def random_shuffling_trial(egrm_results, sriram_results, overlapping_ids):
    # this is for one chromosome, but all overlapping individuals

    # true average overlap
    true_avg_overlap = 0
    for id in overlapping_ids:
        indv_egrm = egrm_results[id]
        indv_sriram = sriram_results[id]
        indv_overlap = calc_overlap(indv_egrm.region_lst, indv_sriram.region_lst)
        sriram_sum = sum([(x.end - x.start) for x in indv_sriram.region_lst])
        true_avg_overlap += indv_overlap/sriram_sum
    true_avg_overlap = true_avg_overlap/len(overlapping_ids)

    # restrict to overlapping for our results, 
    egrm_subset = {k: egrm_results[k] for k in overlapping_ids}

    # shuffle many times and compute avg overlap
    avg_overlap_dist = []
    for t in range(RAND_TRIALS):
        egrm_shuffled = shuffle_indvs(egrm_subset)

        # compute avg overlap
        avg_overlap = 0
        for id in egrm_shuffled:
            indv_egrm = egrm_shuffled[id]
            indv_sriram = sriram_results[id]
            indv_overlap = calc_overlap(indv_egrm.region_lst, indv_sriram.region_lst)
            sriram_sum = sum([(x.end - x.start) for x in indv_sriram.region_lst])
            avg_overlap += indv_overlap/sriram_sum

        avg_overlap = avg_overlap/len(egrm_shuffled)
        avg_overlap_dist.append(avg_overlap)

    # compute p-value
    mean, std = np.mean(avg_overlap_dist), np.std(avg_overlap_dist)
    print("mean, std, true", mean, std, true_avg_overlap)
    pvalue = compute_pvalue(true_avg_overlap, mean, std)
    return pvalue

def make_random_indv(region_lst, chrom_len):
    # make calls for a random indv using the same set of regions
    # note: we don't check that none of the regions overlap...
    rand_region_lst = []
    for region in region_lst:
        region_len = region.end - region.start
        new_start = random.randrange(chrom_len)
        new_end = new_start + region_len
        new_region = Region(region.chrom, new_start, new_end, region.prob)
        rand_region_lst.append(new_region)
    
    return rand_region_lst

def random_trials(proposed_region_lst, target_region_lst, chrom_len):
    # make 1000 random haps to get a distribution (for finding p-value later)
    overlap_dist = []
    for i in range(RAND_TRIALS):
        rand_region_lst = make_random_indv(proposed_region_lst, chrom_len)
        rand_overlap = calc_overlap(rand_region_lst, target_region_lst)
        overlap_dist.append(rand_overlap)
    return overlap_dist

def compute_pvalue(value, mean, std):
    # right now this is one-sided (greater than mean)
    test_stat = (value-mean)/std
    #pvalue, _ = quad(norm.pdf, test_stat, float('inf')) # other way
    pvalue = 1 - norm.cdf(test_stat)
    return pvalue

def liftOver(region_lst):
    # convert from hg19 to hg38 using cumbersome procedure..
    #print("num before", len(region_lst))
    
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

    #print("num after", len(new_regions))
    return new_regions

def read_egrm(target_chrom, pred_filename, id_filename):
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

    # merge consecutive regions
    for id in egrm_results:
        egrm_results[id].merge_consecutive()
    return egrm_results

def read_sriram(target_chrom, pred_filename, id_filename):
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

            if id not in sriram_results:
                sriram_results[id] = Individual(target_chrom)
            
            sriram_results[id].add_region(region)
    
    pred_file.close()

    # convert to hg38 (liftOver)
    for id in sriram_results:
        #print("liftOver", id)
        indv = sriram_results[id]
        indv.region_lst = liftOver(indv.region_lst)

    return sriram_results

def read_ibdmix(target_chrom, pred_filename, target_pop):
    pred_file = open(pred_filename,'r')
    ibdmix_results = {}

    for line in pred_file:
        tokens = line.split()
        chrom = tokens[0][3:]
        pop = tokens[6]
        if chrom == target_chrom and pop == target_pop:
            
            id = tokens[8]
            start = int(tokens[1])
            end = int(tokens[2])
            prob = float(tokens[3]) # LOD not prob...

            region = Region(chrom, start, end, prob)

            if id not in ibdmix_results:
                ibdmix_results[id] = Individual(target_chrom)
            
            ibdmix_results[id].add_region(region)

    pred_file.close()
    return ibdmix_results

def read_skov(target_chrom, pred_filename, target_pop):
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
            # TODO tokens[5] is none, Nea, Den... should we drop none or have threshold 0.8?

            region = Region(chrom, start, end, prob)

            if id not in skov_results:
                skov_results[id] = Individual(target_chrom)
            
            skov_results[id].add_region(region)
    
    pred_file.close()
    return skov_results

def calc_overlap(regions1, regions2):
    # total overlap between all regions of 1 and 2
    total_overlap = 0
    for regionA in regions1:
        for regionB in regions2:
            overlap = overlap_one_region(regionA, regionB)
            if overlap > 0:
                total_overlap += overlap

    return total_overlap

def overlap_one_region(regionA, regionB):
    assert regionA.chrom == regionB.chrom
    a1 = regionA.start
    a2 = regionA.end

    b1 = regionB.start
    b2 = regionB.end

    # calcuate overlap (negative if no overlap)
    return min(a2,b2) - max(a1,b1)

################################################################################
# MAIN
################################################################################

def avg_frac_nea(results_dict):
    """average frac nea over all indvs, not just those overlapping with other methods
    right now this is per chrom"""
    total = 0
    for id in results_dict:
        total += results_dict[id].frac_nea()
    return total/len(results_dict)

def pairwise_overlap(resultsA, resultsB):
    """from two methods, each with a results dictionary (key: indv ID, value: indv)
    compute their overlap. Right now this is per chrom"""

    overlapping_ids = set(resultsA.keys()).intersection(set(resultsB.keys()))
    num_indv = len(overlapping_ids)
    print("num overlapping", num_indv)

    #for id in resultsA:
    #    if id not in resultsB:
    #        print("in our data but not 1000g?", id)
    #input('enter')

    #chr_len = CHROM_DICT[CHR]

    # meta stats
    avg_overlap = 0
    frac_pvalue_sig = 0

    # go through each individual where we overlap
    for id in overlapping_ids:
        #print("\nstarting hap", id)
        indv_egrm = resultsA[id]
        indv_sriram = resultsB[id]
        #print("egrm", indv_egrm)
        #print("sriram", indv_sriram)
        indv_overlap = calc_overlap(indv_egrm.region_lst, indv_sriram.region_lst)

        sriram_sum = sum([(x.end - x.start) for x in indv_sriram.region_lst])

        frac_overlap = indv_overlap/sriram_sum
        #print("overlap/sriram indv", frac_overlap)
        avg_overlap += frac_overlap

        # random testing per individual
        '''overlap_dist = random_trials(indv_egrm.region_lst, indv_sriram.region_lst, chr_len)
        rescaled_dist = [x/sriram_sum for x in overlap_dist]
        mean, std = np.mean(rescaled_dist), np.std(rescaled_dist)
        print("mean, std", mean, std)
        pvalue = compute_pvalue(frac_overlap, mean, std)
        if pvalue <= SIG_THRESH:
            frac_pvalue_sig += 1
        print("pvalue", pvalue)'''

    avg_overlap = avg_overlap/num_indv
    print("avg indv overlap fraction", avg_overlap)
    #shuffle_pvalue = random_shuffling_trial(resultsA, resultsB, overlapping_ids)
    #print("shuffle indv pvalue", shuffle_pvalue)
    #print("frac sig p-values", frac_pvalue_sig/num_indv)

def one_chrom(CHR, egrm_pred_filename, egrm_id_filename, sriram_pred_filename, sriram_id_filename, ibdmix_pred_filename):
    # dictionaries of individual (key) : regions (value)
    egrm_results = read_egrm(CHR, egrm_pred_filename, egrm_id_filename)
    sriram_results = read_sriram(CHR, sriram_pred_filename, sriram_id_filename)
    ibdmix_results = read_ibdmix(CHR, ibdmix_pred_filename, POP)

    print("\nstarting chrom", chr)
    print("egrm num indvs", len(egrm_results), "frac nea", avg_frac_nea(egrm_results))
    print("sriram num indvs", len(sriram_results), "frac nea", avg_frac_nea(sriram_results))
    print("ibdmix num indvs", len(ibdmix_results), "frac nea", avg_frac_nea(ibdmix_results))

    # compare all pairs! TODO make matrix
    pairwise_overlap(egrm_results, sriram_results)
    pairwise_overlap(egrm_results, ibdmix_results)
    pairwise_overlap(ibdmix_results, sriram_results)

if __name__ == "__main__":
    for chr_int in range(1,23):
        chr = str(chr_int)
        egrm_pred_filename = egrm_pred_path + POP + "_chr" + chr + ".pred"
        sriram_pred_filename = sriram_pred_path + "chr-" + chr + ".thresh-90.length-0.00.haplotypes.gz"
        one_chrom(str(chr_int), egrm_pred_filename, egrm_id_filename, sriram_pred_filename, sriram_id_filename, ibdmix_pred_filename)
