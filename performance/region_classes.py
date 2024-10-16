"""
Helper classes (Individual and Region) for comparing methods
Author: Sara Mathieson
Date: 10/16/24
"""

import numpy as np

################################################################################
# HELPERS
################################################################################

def read_chrom_lengths():
    arr = np.loadtxt("../gnomad/hg38_chrom_lengths.tsv", dtype='int', delimiter="\t", skiprows=1)
    chrom_dict = {}
    for chr in range(1,23):
        assert arr[chr-1][0] == chr
        chrom_dict[str(chr)] = arr[chr-1][1]
    return chrom_dict

CHROM_DICT = read_chrom_lengths()

################################################################################
# CLASSES
################################################################################

class Individual:

    def __init__(self, chrom):
        # add population, id number?
        self.chrom = chrom
        self.region_lst = []

    def add_region(self, region):
        if region not in self.region_lst:
            self.region_lst.append(region)
        else:
            #print("region already in region_lst!!")
            #input('enter')
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

        #print("merged", len(merged_idx_lst), len(self.region_lst))
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
    
    def collapse_haps(self):
        # if homozygous, collapse
        #print("starting to collapse!")
        #print("num regions", len(self.region_lst))
        num_regions = len(self.region_lst)

        # find overlapping indices
        overlapping_indices = []
        all_indices = set()
        for i in range(num_regions):
            for j in range(i+1, num_regions):
                regionA = self.region_lst[i]
                regionB = self.region_lst[j]
                if regionA.overlap(regionB) > 0:
                    #print("self overlap!!", i, j, regionA, regionB)
                    all_indices.add(i)
                    all_indices.add(j)
                    found = False
                    for group in overlapping_indices:
                        if i in group or j in group:
                            group.add(i)
                            group.add(j)
                            found = True
                    if not found:
                        overlapping_indices.append(set([i,j]))
        #print(overlapping_indices)

        merged_regions = []
        for group in overlapping_indices:
            # make new region with min start and max end
            start = min([self.region_lst[i].start for i in group])
            end = max([self.region_lst[i].end for i in group])
            chrom = self.region_lst[list(group)[0]].chrom # should all be the same
            prob = np.mean([self.region_lst[i].prob for i in group]) # avg
            new_region = Region(chrom, start, end, prob)
            merged_regions.append(new_region)
        
        # add back in all the rest
        for i in range(num_regions):
            if i not in all_indices:
                merged_regions.append(self.region_lst[i])

        self.region_lst = merged_regions
        #print("after", len(self.region_lst))

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

    def overlap(self, other):
        if self.chrom != other.chrom:
            print(self.chrom, other.chrom)
        assert self.chrom == other.chrom
        a1 = self.start
        a2 = self.end

        b1 = other.start
        b2 = other.end

        # calcuate overlap (negative if no overlap)
        return min(a2,b2) - max(a1,b1)

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