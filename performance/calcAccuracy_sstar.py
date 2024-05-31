import pandas as pd
import sys
import pybedtools
import pickle
import numpy as np
from tqdm import tqdm

model_type = sys.argv[1]
seed = sys.argv[2]

x = []
y = []
true_positives_list = []
total_inferred_tracts_list = []
total_true_tracts_list = []

def cal_accuracy(true_tracts, inferred_tracts):
    """
    Description:
        Helper function for calculating accuracy.

    Arguments:
        true_tracts str: Name of the BED file containing true introgresssed tracts.
        inferred_tracts str: Name of the BED file containing inferred introgressed tracts.

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """
    truth_tracts = pybedtools.bedtool.BedTool(true_tracts).sort().merge().intersect(pybedtools.BedTool("interval_25k_175k.bed"))
#     print(truth_tracts)
    inferred_tracts =  pybedtools.bedtool.BedTool(inferred_tracts).merge()

    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

    return true_positives, total_inferred_tracts, total_true_tracts


# iterate cut off
for cutoff in tqdm([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999,0.9999]):
# for cutoff in tqdm([0.1,0.2,0.3]):
    true_positives = 0
    total_inferred_tracts = 0
    total_true_tracts = 0
    # iterate runs
    for run in range(100):
        # iterate individuals
        df = pd.read_csv(f'/scratch1/jcahoon/argml/data/distribution_{model_type}_validate/sstar_score/{seed}/{cutoff}/200kb_{seed}_{run}_{cutoff}.tracts.bed', delimiter='\t', names=["", "left", "right", "node_sample"])
        for indv in range(56):
            temp = df[df["node_sample"] == f'tsk_{indv}']
            segments=[]
            if temp.empty:
                segments.append([0,0])
            else:
                for ind in temp.index:
                    segments.append([temp['left'][ind], temp['right'][ind]])
            with open(f'/dev/shm/200kb_{seed}_{run}_{indv}.bed', 'w') as o:
                for seg in segments:
                    o.write(f'1\t{int(seg[0])}\t{int(seg[1])}\n')
            tp, tit, ttt = cal_accuracy(f'/scratch1/jcahoon/argml/data/distribution_{model_type}_validate/segments_indv/{seed}/200kb_{seed}_{run}_{indv}.bed',f'/dev/shm/200kb_{seed}_{run}_{indv}.bed')
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
with open(f'pr_pickle/200kb_0_precision_recall_sstar_{model_type}_{seed}.pickle', "wb") as f:
    pickle.dump([x,y, true_positives_list, total_inferred_tracts_list, total_true_tracts_list],f)
