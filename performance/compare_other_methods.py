"""
Compare introgression results from four different methods.
Author: Sara Mathieson
Date: 10/16/24
"""

# python imports
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# our imports
import parse_egrm, parse_sriram, parse_ibdmix, parse_skovhmm

################################################################################
# GLOBALS
################################################################################

COMPARE_PATH = "/homes/smathieson/Documents/arg-ml/"
POP = "CEU"
#RAND_TRIALS = 1000
#SIG_THRESH = 0.05 # significance threshold (p-value)

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

def calc_overlap(regions1, regions2):
    # total overlap between all regions of 1 and 2
    total_overlap = 0
    for regionA in regions1:
        for regionB in regions2:
            overlap = regionA.overlap(regionB)
            if overlap > 0:
                total_overlap += overlap

    return total_overlap

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
    return avg_overlap

def one_chrom(names, all_results):
    for result in all_results:
        print("num indvs", len(result), "frac nea", avg_frac_nea(result))

    # compare all pairs
    num_methods = len(all_results)
    pairwise_overlaps = np.zeros((num_methods, num_methods))
    for i in range(num_methods):
        for j in range(num_methods):
            print(names[i], names[j])
            avg_overlap = pairwise_overlap(all_results[i], all_results[j])
            pairwise_overlaps[i][j] = avg_overlap

    sns.heatmap(pairwise_overlaps, annot=True, xticklabels=names, yticklabels=names, cmap="rocket_r", vmin=0, vmax=1)
    plt.title(POP + " chr" + CHR + " method comparison")
    plt.savefig(COMPARE_PATH + "output/" + POP + "/" + POP + "_chr" + CHR + "_overlaps2.pdf")
    plt.clf()

################################################################################
# MAIN
################################################################################

if __name__ == "__main__":
    # method names
    names = ["EGRM", "Sriram", "IBDmix", "Skov_HMM"]

    for chr_int in range(1,23):
        CHR = str(chr_int)
        print("\nstarting chrom", CHR)
        egrm_results = parse_egrm.egrm_one_chrom(CHR)
        sriram_results = parse_sriram.sriram_one_chrom(CHR)
        ibdmix_results = parse_ibdmix.ibdmix_one_chrom(CHR)
        skovhmm_results = parse_skovhmm.skovhmm_one_chrom(CHR)

        all_results = [egrm_results, sriram_results, ibdmix_results, skovhmm_results]
        one_chrom(names, all_results)