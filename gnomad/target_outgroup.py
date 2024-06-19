"""
Description: write new samples files with both target and outgroup
Usage: python3 target_outgroup.py
Author: Sara Mathieson
Date: 6/19/24
"""

# python imports
import os

################################################################################
# GLOBALS
################################################################################

HEADER = "gnomad_subpops/header.txt"
OUTGROUP = "gnomad_subpops/yri.txt"
N = 56 # num of individuals in the target, num individuals in the outgroup

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
            
def main():

    pop_file_lst = os.listdir("gnomad_subpops")
    pop_file_lst.sort()

    # header and YRI
    header_str = open(HEADER, "r").read()
    yri_str = get_yri(header_str)

    # go through each population
    for pop_file in pop_file_lst:
        
        pop = pop_file[:3].upper()
        if pop != "YRI":
            old_pop_filename = "gnomad_subpops/" + pop_file
            new_pop_filename = "gnomad_subpops/" + pop_file[:-4] + "_yri.txt"
            print(old_pop_filename, new_pop_filename)

            old_pop_file = open(old_pop_filename, "r")
            new_pop_file = open(new_pop_filename, "w")
            count = 0
            for line in old_pop_file:
                id = line.strip()
                if id in header_str and count < N:
                    new_pop_file.write(id + "\n")
                    count += 1

            # close files
            old_pop_file.close()
            
            if count != N:
                print("ERROR: not enough individuals in " + pop)
        
        new_pop_file.write(yri_str)
        new_pop_file.close()
        input('enter')


if __name__ == "__main__":
    main()