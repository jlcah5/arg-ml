"""
Description: write new samples files with both target (1st) and outgroup (2nd)
Usage: python3 target_outgroup.py
Author: Sara Mathieson
Date: 6/19/24
"""

# python imports
import os
import sys

################################################################################
# GLOBALS
################################################################################

HEADER = "gnomad_subpops/header.txt"
OUTGROUP = "gnomad_subpops/yri.txt"
#SRIRAM_IDS = sys.argv[1]
N = 56 # num of individuals in the target, num individuals in the outgroup
#TODO stopped here: make ids that match sriram and overlap slightly

################################################################################
# MAIN
################################################################################

def get_yri(header_str):
    yri_filename = "gnomad_subpops/yri.txt"
    yri_file = open(yri_filename, "r")
            
    count = 0
    yri_str = ""
    for line in yri_file:
        id = line.strip()
        if id in header_str and count < N:
            yri_str += id + "\n"
            count += 1

    # close files
    yri_file.close()
    return yri_str  

'''def read_sriram_ids(sriram_filename):
    sriram_file = open(sriram_filename,'r')
    sriram_ids = set()
    for line in sriram_file:
        tokens = line.split()
        sriram_ids.add(tokens[0])
    return list(sriram_ids)'''

def main():

    pop_file_lst = os.listdir("gnomad_subpops")
    pop_file_lst.sort()

    # header and YRI
    header_str = open(HEADER, "r").read()
    yri_str = get_yri(header_str)

    # go through each population
    for pop_file in pop_file_lst:
        
        pop = pop_file[:3].upper()
        #if pop != "YRI" and pop != "HEA": # outgroup or header.txt
        if pop == "CEU" and "yri" not in pop_file: # just do for ceu right now
            old_pop_filename = "gnomad_subpops/" + pop_file
            
            new_pop_filename = "gnomad_subpops/" + pop_file[:-4] + "_yri2.txt"
            print(old_pop_filename, new_pop_filename)

            old_pop_file = open(old_pop_filename, "r")
            new_pop_file = open(new_pop_filename, "w")
            #count = 0
            valid_ids = []
            for line in old_pop_file:
                id = line.strip()
                if id in header_str: # I think this check is unnecessary but just in case
                    valid_ids.append(id)

            # close files
            old_pop_file.close()

            print("num valid", len(valid_ids))
            for i in range(N,2*N):
                new_pop_file.write(valid_ids[i] + "\n")
            
            #if count != N:
            if len(valid_ids) < 2*N:
                print("ERROR: not enough individuals in " + pop)
        
            new_pop_file.write(yri_str) # YRI is second
            new_pop_file.close()


if __name__ == "__main__":
    main()