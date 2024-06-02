import pandas as pd
import tskit
import numpy as np
import tqdm
import pickle
import pybedtools
import sys
from utils import calc_accuracy

# EXAMPLE USAGE: python calcAccuracy-nn.py growth 4 200 transformer 600000
# TO-DO: edit cutoff if altering cutoff in calcTractThresh.sh
# USER INPUT
model_type = sys.argv[1] # growth or nonGrowth
seed = sys.argv[2] # seed used for simulation
len_prefix = sys.arv[3] # length in kb
method = sys.argv[4] # transformer or cnn
train_size =  int(int(sys.argv[5])/50000) # number of training samples, in multiples of 50k

# custum variables
cutoff_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.9,0.95, 0.99,0.995, 0.999] # list of tresholds, customizable
numInstances = 100
num_haplotypes = 56

# open prediction outputs and true labels
with open(f'~/arg-ml/ml/pr_curve/Distribution_{method}_50k_sw_{model_type}_{iter_num}_{train_size*50}K_y_out.npy', "rb") as f:
    pred_arr= np.array(pickle.load(f))
with open(f'~/arg-ml/ml/pr_curve/Distribution_{method}_50k_sw_{model_type}_{iter_num}_{train_size*50}K_y_test.npy', "rb") as f:
    test_arr= pickle.load(f)
    
x=[]
y=[]
true_positives_list = []
total_inferred_tracts_list = []
total_true_tracts_list = []

for cutoff in cutoff_list:
    total_inferred_tracts = 0
    total_true_tracts = 0
    true_positives = 0 
    for tree in tqdm.tqdm(range(num_instances)):
        df = pd.read_csv(f'~/arg-ml/data/distribution_{model_type}/segments/{iter_num}/{len_prefix}kb_{seed}_{tree}.csv')
        for hap in range(num_haplotypes): # iterate the haplotypes
            # get true seg
            true=[]
            df_subset = df[(df['node_sample']==2*hap) | (df['node_sample']==(2*hap+1))]
            if not df_subset.empty:
                for row in df_subset.iterrows():
                    row=row[1]
                    true.append((int(row.left), int(row.right)))
            start_hap = hap*2
            end_hap = hap*2+2

            start_tree = tree*16
            end_tree = tree*16+16
            pred_sw = pred_arr[start_tree:end_tree, start_hap:end_hap]
            pred = []
            in_seg = False
            start_pos = 0 
            end_pos = 0
            for pos in range(16):
                hit = pred_sw[pos,:][0] > cutoff or pred_sw[pos,:][1] > cutoff
                if hit and not in_seg:
                    in_seg = True
                    start = int((pos*10000+ pos*10000 + 50000)/2)
                elif not hit and in_seg:
                    in_seg = False
                    end = int((pos*10000+ pos*10000 + 50000)/2)
                    pred.append((start,end))
                    start = 0
                    end = 0
                elif hit and pos == 15:
                    end = int((pos*10000+ pos*10000 + 50000)/2)
                    pred.append((start,end))
            true_file = f'/dev/shm/{len_prefix}kb_0_{tree}_{method}_{model_type}_{hap}_{seed}_true.bed'
            pred_file = f'/dev/shm/{len_prefix}kb_0_{tree}_{method}_{model_type}_{hap}_{seed}_pred.bed'
            true.sort(key=lambda x:x[0])
            pred.sort(key=lambda x:x[0])
            
            # write to bed for accuracy calculation
            with open(true_file, 'w') as o:
                for t in true:
                    o.write(f'1\t{t[0]}\t{t[1]}\n')
            with open(pred_file, 'w') as o:
                for t in pred:
                    o.write(f'1\t{t[0]}\t{t[1]}\n')
            it, tt, tp = cal_accuracy(true_file,pred_file)
            total_inferred_tracts += it
            total_true_tracts += tt

    if float(total_inferred_tracts) == 0: precision = np.nan
    else: precision = true_positives / float(total_inferred_tracts) * 100
    if float(total_true_tracts) == 0: recall = np.nan
    else: recall = true_positives / float(total_true_tracts) * 100
    x.append(precision)
    y.append(recall)
    true_positives_list.append(true_positives)
    total_inferred_tracts_list.append(total_inferred_tracts)
    total_true_tracts_list.append(total_true_tracts)
    if recall == 0:
        break
print(x,y)

# save output
with open(f'pr_pickle/{len_prefix}kb_precision_recall_Distribution_{method}_{model_type}_{seed}_{train_size*50}K.pickle', "wb") as f:
    pickle.dump([x,y, true_positives_list, total_inferred_tracts_list, total_true_tracts_list],f)