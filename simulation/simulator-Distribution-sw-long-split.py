import sys
import msprime
import numpy as np
from numpy.random import default_rng
import pandas as pd
import tskit
from tqdm import tqdm
import subprocess
import demes
from utils import get_introgressed_tracts

# input params
base_seed = int(sys.argv[1])
num = int(sys.argv[2])
prefix_length  = int(sys.argv[3])
model_number = sys.argv[4]
model_type = sys.argv[5]
dir_name = sys.argv[6]
window_size = int(sys.argv[7])
step = int(sys.argv[8])


rng=default_rng(base_seed)
n=56
length = prefix_length*1000

source_id = 2
target_id = 3

for i in tqdm(range(num)):

    ts = tskit.load(f'{dir_name}/trees/{base_seed}/{prefix_length}kb_{base_seed}_{i}.tree')

    for start in range(0,length-window_size+step,step):
            df= []
            ts_new = ts.keep_intervals(np.array([[start, start+window_size]]), simplify=False).trim()
            if len(get_introgressed_tracts(ts_new, source_id, target_id)) > 0:
                for m in ts_new.migrations():
                    introgressed_set=set()
                    if m.dest == source_id and m.source == target_id:
        #                 print(m)
                        temp=[]
                        for tree in ts_new.trees():
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
                if len(df[['node_nea','left', 'right']].value_counts().keys()) != len(get_introgressed_tracts(ts_new, source_id, target_id)): # double check all segments are accounted for
                    print("-----")
                    print("WRONG")
                    print(df)
                    print(df[['node_nea','left', 'right']].value_counts().keys())
                    print(get_introgressed_tracts(ts_new, source_id, target_id))
                    print("-----")
            else:
                df=pd.DataFrame([], columns=['node_sample', 'node_nea', 'left', 'right'])
        #         print(df)
            # with write to vcf
            with open('{}/trees_vcf/{}/{}kb_{}_{}_{}.infer.vcf'.format(dir_name, base_seed, prefix_length, base_seed,i, start), "w") as f:
                ts_new.write_vcf(output=f)

            # drump tree file
            ts_new.dump('{}/trees/{}/{}kb_{}_{}_{}.infer.tree'.format(dir_name, base_seed, prefix_length,base_seed,i,start)) 

            # output df
            df.to_csv('{}/segments/{}/{}kb_{}_{}_{}.infer.csv'.format(dir_name, base_seed, prefix_length,base_seed,i,start), index=False)
            #convert to egrm
            tree_name=f'{dir_name}/trees/{base_seed}/{prefix_length}kb_{base_seed}_{i}_{start}.infer.tree'
            egrm_name=f'{dir_name}/egrm/{base_seed}/{prefix_length}kb_{base_seed}_{i}_{start}.infer.input'
            subprocess.run([f'trees2egrm', '--output-format', 'numpy', f'{tree_name}', '--c', '--haploid', '--output', f'{egrm_name}'],
                          stdout = subprocess.DEVNULL,
                           stderr = subprocess.DEVNULL) 

            # cleaning log files
            subprocess.run(['rm', f'{dir_name}/egrm/{base_seed}/{prefix_length}kb_{base_seed}_{i}_{start}.infer.input_mu.npy', f'{dir_name}/egrm/{base_seed}/{prefix_length}kb_{base_seed}_{i}_{start}.infer.input.log.log']) 

            label = np.zeros(n*4) # 2*(2n)
            if len(df) > 0:
                label[list(df['node_sample'].unique())] = 1
            np.save('{}/egrm/{}/{}kb_{}_{}_{}.infer.output'.format(dir_name,base_seed,prefix_length,base_seed,i,start), label)

