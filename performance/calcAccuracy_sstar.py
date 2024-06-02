import pandas as pd
import sys
import pybedtools
import pickle
import numpy as np
from tqdm import tqdm
from utils import calc_accuracy

# EXAMPLE USAGE: python calcAccuracy-sstar.py growth 4 200
# TO-DO: edit cutoff if altering cutoff in calcTractThresh.sh
# USER INPUT
model_type = sys.argv[1] # growth or nonGrowth
seed = sys.argv[2] # seed used for simulation
len_prefix = sys.arv[3] # length in kb

# list to track output
x = []
y = []
true_positives_list = []
total_inferred_tracts_list = []
total_true_tracts_list = []

# iterate cut off
for cutoff in tqdm([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999,0.9999]):
    true_positives = 0
    total_inferred_tracts = 0
    total_true_tracts = 0
    # iterate runs
    for run in range(100):
        # iterate individuals
        df = pd.read_csv(f'~/arg-ml/data/distribution_{model_type}/sstar_score/{seed}/{cutoff}/{len_prefix}kb_{seed}_{run}_{cutoff}.tracts.bed', delimiter='\t', names=["", "left", "right", "node_sample"])
        for indv in range(56):
            temp = df[df["node_sample"] == f'tsk_{indv}']
            segments=[]
            if temp.empty:
                segments.append([0,0])
            else:
                for ind in temp.index:
                    segments.append([temp['left'][ind], temp['right'][ind]])
            with open(f'/dev/shm/{len_prefix}kb_{seed}_{run}_{indv}.bed', 'w') as o:
                for seg in segments:
                    o.write(f'1\t{int(seg[0])}\t{int(seg[1])}\n')
            tp, tit, ttt = cal_accuracy(f'/arg-ml/data/distribution_{model_type}/segments_indv/{seed}/{len_prefix}kb_{seed}_{run}_{indv}.bed',f'/dev/shm/{len_prefix}kb_{seed}_{run}_{indv}.bed')
            true_positives += tp
            total_inferred_tracts += tit
            total_true_tracts += ttt
    # calculate precision recall
    if float(total_inferred_tracts) == 0: precision = np.nan
    else: precision = true_positives / float(total_inferred_tracts) * 100
    if float(total_true_tracts) == 0: recall = np.nan
    else: recall = true_positives / float(total_true_tracts) * 100
    x.append(precision)
    y.append(recall)
    true_positives_list.append(true_positives)
    total_inferred_tracts_list.append(total_inferred_tracts)
    total_true_tracts_list.append(total_true_tracts)
print(x,y)

# save output
with open(f'pr_pickle/{len_prefix}kb_precision_recall_sstar_{model_type}_{seed}.pickle', "wb") as f:
    pickle.dump([x,y, true_positives_list, total_inferred_tracts_list, total_true_tracts_list],f)
