import sys
import msprime
import numpy as np
from numpy.random import default_rng
import pandas as pd
import tskit
from tqdm import tqdm
import subprocess
import demes

def get_introgressed_tracts(ts, src_name, tgt_name):
    """
    Description:
        Outputs true introgressed tracts from a tree-sequence into a BED file.

    Arguments:
        ts tskit.TreeSequence: Tree-sequence containing introgressed tracts.
        chr_name int: Name of the chromosome.
        src_name str: Name of the source population.
        tgt_name str: Name of the target population.
        output string: Name of the output file.
    """
    introgressed_tracts = []
    for m in ts.migrations():
#         print(m)
        if m.dest == source_id and m.source == target_id: introgressed_tracts.append((int(m.left), int(m.right)))

    introgressed_tracts.sort(key=lambda x:x[0])
    return  introgressed_tracts

# input params
base_seed = int(sys.argv[1])
num = int(sys.argv[2])
prefix_length  = int(sys.argv[3])
model_number = sys.argv[4]
model_type = sys.argv[5]
dir_name= sys.argv[6] # output directory


rng=default_rng(base_seed)
n=56 # number of individuals
length = prefix_length*1000

# setting up directory
for i in ['trees','trees_vcf','egrm','segments']:
    subprocess.run(['mkdir', f'{dir_name}/{i}/{base_seed}'],
                   stdout = subprocess.DEVNULL,
                   stderr = subprocess.DEVNULL)  

# loading demographic model
graph = demes.load(f'demographic_models/distribution_{model_type}/yaml/model_{model_number}.yaml')
demography = msprime.Demography.from_demes(graph)
demography.sort_events()
source_id = 2
target_id = 3

for i in tqdm(range(num)):
    seed = rng.integers(1,high=2**32) 
    ts = msprime.sim_ancestry(samples={'ANC':0, "SRC": 0 , "TAR": n, "OUT": n}, random_seed=seed, demography=demography, sequence_length=length, ploidy=2, recombination_rate = 1.25e-8, record_migrations=True)
    ts= msprime.sim_mutations(ts, rate = 1.29e-8, discrete_genome = False, random_seed=seed)
    df= []
    if len(get_introgressed_tracts(ts, source_id, target_id)) > 0:
        for m in ts.migrations():
            introgressed_set=set()
            if m.dest == source_id and m.source == target_id:
#                 print(m)
                temp=[]
                for tree in ts.trees():
                    if tree.interval.left >= m.left and tree.interval.right <= m.right:
                        if tree.is_sample(m.node): # check if migration happens to sample node, then add to introgressed list
                            introgressed_set.add(m.node)
                        else: # BFS to find sampled nodes
                            child_list=tree.children(m.node)
                            while len(child_list) > 0:
                                curr = child_list[0]
                                if tree.is_sample(curr):
                                    introgressed_set.add(curr)
                                child_list = child_list[1:]
                                child_list += tree.children(curr)
#                 print(introgressed_set)
                for indv in introgressed_set: # create df rows
                    temp.append([int(indv), m.node, m.left, m.right])
                df.append(pd.DataFrame(temp, columns=['node_sample', 'node_nea', 'left', 'right']))
    if len(df) > 0:
        df = pd.concat(df)
        if len(df[['node_nea','left', 'right']].value_counts().keys()) != len(get_introgressed_tracts(ts, source_id, target_id)): # double check all segments are accounted for
            print("-----")
            print("WRONG")
            print(df)
            print(df[['node_nea','left', 'right']].value_counts().keys())
            print(get_introgressed_tracts(ts, source_id, target_id))
            print("-----")
    else:
        df=pd.DataFrame([], columns=['node_sample', 'node_nea', 'left', 'right'])
#         print(df)
    # with write to vcf
    with open('{}/trees_vcf/{}/{}kb_{}_{}.vcf'.format(dir_name, base_seed, prefix_length, base_seed,i), "w") as f:
        ts.write_vcf(output=f)
        
    # drump tree file
    ts.dump('{}/trees/{}/{}kb_{}_{}.tree'.format(dir_name, base_seed, prefix_length,base_seed,i)) 
    
    # output df
    df.to_csv('{}/segments/{}/{}kb_{}_{}.csv'.format(dir_name, base_seed, prefix_length,base_seed,i), index=False)
    

