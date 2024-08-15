"""
Evaluate predictions for real data (right now through B maps)
https://github.com/sellalab/HumanLinkedSelectionMaps/tree/master/Bmaps
phastCons_bestfit.tar.gz
Idea: map from 10kb regions to avg nea prob and b-stat values, then bin the
bstats (quintiles) and average preds within each one
Author: Sara Mathieson
Date: 7/23/24
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

PRED_FOLDER = sys.argv[1]
BMAP_FOLDER = sys.argv[2]
OUT_FOLDER  = sys.argv[3]

POPS = ["CEU", "GBR"]
REGION_LEN = 10000 # split 50kb inference into 10kb regions
NUM_HAPS = 112

class Region:

    def __init__(self, chr, start, end):
        self.chr = chr # int
        self.start = start
        self.end = end
        self.preds = []

    def find_bstat(self, bstat_lst):
        # bstat list should correspond to relevant chromosome
        if self.chr != bstat_lst[0][3]:
            print(self.chr, bstat_lst[0][3])
        assert self.chr == bstat_lst[0][3] # check chrom
        region_start_idx, start_inside = binary_search(self.start, bstat_lst)
        region_end_idx, end_inside = binary_search(self.end, bstat_lst)

        region_start = bstat_lst[region_start_idx]
        region_end = bstat_lst[region_end_idx]

        # weighted average of bstats
        if (not start_inside) or (not end_inside):
            #print("Error: region not in bmap")
            pass
        else:
            bstats = []
            lengths = []
            for i in range(region_start_idx, region_end_idx+1):
                bstat = bstat_lst[i][2]
                # i = region_start_idx
                if i == region_start_idx:
                    length = min(self.end, bstat_lst[i][1]) - self.start

                # region_start_idx < i < region_end_idx
                elif region_start_idx < i < region_end_idx:
                    length = bstat_lst[i][1] - bstat_lst[i][0]

                # i = region_end_idx
                elif i == region_end_idx:
                    length = self.end - max(self.start, bstat_lst[i][0])

                else:
                    print("Error!")
                    input('enter')

                bstats.append(bstat)
                lengths.append(length)

            # adjut for end exclusive
            for i in range(region_end_idx-region_start_idx):
                lengths[i] += 1
            #if sum(lengths) != REGION_LEN:
            #    print("Not really error but something going on with lengths..")
            #    print(sum(lengths))
            return np.average(bstats, weights=lengths)

    def avg_preds(self):
        if len(self.preds) != NUM_HAPS*5:
            #print("Error: region on end of chrom")
            pass
        else:
            # take middle 3 preds, averaged over all indvs
            retained_preds = self.preds[NUM_HAPS:-NUM_HAPS]
            assert len(retained_preds) == NUM_HAPS*3
            return sum(retained_preds)/len(retained_preds)

def read_bmap(filename):
    """Read bstats from bed file"""

    bstat_lst = []
    f = open(filename,'r')

    for line in f:
        tokens = line.split()
        chrom = int(tokens[0][3:].split("_")[0])
        start = int(tokens[1])
        end = int(tokens[2]) # inclusive
        bstat = float(tokens[3])/1000 # rescale

        bstat_lst.append([start,end,bstat,chrom])

    f.close()
    return bstat_lst

def binary_search(q, lst):
    low = 0
    high = len(lst)-1

    while low <= high:

        mid = (low+high)//2
        if lst[mid][0] <= q <= lst[mid][1]: # inside region
            return mid, True
        elif q < lst[mid][0]:
            high = mid-1
        else:
            low = mid+1

    return mid, False # something close

def bmap2bed(chr):
    # convert bmap to bed format
    in_filename = BMAP_FOLDER + "/chr" + str(chr) + ".bmap.txt"
    out_filename = BMAP_FOLDER + "_bed/chr" + str(chr) + ".bmap.bed"
    print(in_filename, out_filename)
    in_file = open(in_filename, 'r')
    out_file = open(out_filename, 'w')

    start = 0
    for line in in_file:
        tokens = line.split()
        bstat = tokens[0]
        end = start + int(tokens[1]) - 1
        new_line = " ".join(["chr" + str(chr), str(start), str(end), bstat])
        out_file.write(new_line + "\n")

        # update start
        start = end + 1
    
    in_file.close()
    out_file.close()

def read_pred_file(filename):
    pred_file = open(filename)
    region_dict = {} # key is string i.e. "chr3:10000-20000"
    for line in pred_file:
        tokens = line.split()
        chr = int(tokens[0])
        start = int(tokens[1])
        end = int(tokens[2])
        pred = float(tokens[3])

        # key list (split into 10kb regions)
        for begin in range(start, end, REGION_LEN):
            # add to dictionary
            key = tokens[0] + ":" + str(begin) + "-" + str(begin+REGION_LEN)
            if key not in region_dict:
                region = Region(chr, begin, begin+REGION_LEN)
                region_dict[key] = region

            # add to preds regardless
            region_dict[key].preds.append(pred)

    return region_dict

def get_all_preds_bstats(pred_file, bmap_lst):
    # dictionary for all regions
    region_dict = read_pred_file(pred_file)

    all_bstats = []
    all_preds = []
    for key in region_dict:
        region = region_dict[key]
        bstat = region.find_bstat(bmap_lst)
        avg_pred = region.avg_preds()
        #print("bstat", bstat, "avg pred", avg_pred)
        if bstat is not None and avg_pred is not None:
            all_bstats.append(bstat)
            all_preds.append(avg_pred)

    return all_bstats, all_preds

def quintile_preds(genome_bstats, genome_preds):
    # visualization
    assert len(genome_bstats) == len(genome_preds)
    quintiles = np.quantile(genome_bstats, [0.2, 0.4, 0.6, 0.8])
    print("bstat quintiles", quintiles)
    quintile_preds = [[],[],[],[],[]]
    for i in range(len(genome_preds)):
        bstat = genome_bstats[i]
        pred = genome_preds[i]

        # 5 cases TODO make cleaner
        if 0 <= bstat < quintiles[0]:
            quintile_preds[0].append(pred)
        elif quintiles[0] <= bstat < quintiles[1]:
            quintile_preds[1].append(pred)
        elif quintiles[1] <= bstat < quintiles[2]:
            quintile_preds[2].append(pred)
        elif quintiles[2] <= bstat < quintiles[3]:
            quintile_preds[3].append(pred)
        elif quintiles[3] <= bstat <= 1:
            quintile_preds[4].append(pred)
    
    print("num in each bin", [len(x) for x in quintile_preds])
    return [np.mean(x) for x in quintile_preds]

if __name__ == "__main__":
    # this was for converting bmap to bed (already done)
    #for chr in range(1,23):
    #    bmap2bed(chr)

    for POP in POPS:
        print("\nstarting", POP)

        genome_bstats = []
        genome_preds = []
        for CHR in range(1,23):

            # read bstat file
            bmap_lst = read_bmap(BMAP_FOLDER + "chr" + str(CHR) + ".bmap_hg38.bed")

            # add each chrom
            pred_file = POP + "_chr" + str(CHR) + ".pred"
            print(PRED_FOLDER + POP + "/" + pred_file)
            chr_bstats, chr_preds = get_all_preds_bstats(PRED_FOLDER + POP + "/" + pred_file, bmap_lst)
            genome_bstats.extend(chr_bstats)
            genome_preds.extend(chr_preds)

        binned_preds = quintile_preds(genome_bstats, genome_preds)
        print(binned_preds)
        plt.plot(range(1,6), [np.mean(x) for x in binned_preds], '-o', label=POP)
    
    plt.xlabel("B statistic (binned)")
    plt.ylabel("avg pred introgression")
    plt.legend()
    plt.savefig(OUT_FOLDER + "bstat_" + "_".join(POPS) + ".pdf")
