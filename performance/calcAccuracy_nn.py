import pandas as pd
import tskit
import numpy as np
import tqdm
import pickle
import pybedtools
import sys

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
    inferred_tracts =  pybedtools.bedtool.BedTool(inferred_tracts).sort().merge()

    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

    return  total_inferred_tracts, total_true_tracts, true_positives

model = sys.argv[1]
method = sys.argv[2]
model_type = sys.argv[3]
iter_num = sys.argv[4]
train_size =  int(int(sys.argv[5])/50000)
print(model, method, model_type, iter_num, train_size)
with open(f'/project/lgarcia1_998/jcahoon/argml/ml/pr_curve/{model}_{method}_50k_sw_{model_type}_{iter_num}_{train_size*50}K_y_out.npy', "rb") as f:
    pred_arr= np.array(pickle.load(f))
with open(f'/project/lgarcia1_998/jcahoon/argml/ml/pr_curve/{model}_{method}_50k_sw_{model_type}_{iter_num}_{train_size*50}K_y_test.npy', "rb") as f:
    test_arr= pickle.load(f)
# get calculate accuracy

# print(pred_arr[0])
# print(pred_arr.shape)
x=[]
y=[]
true_positives_list = []
total_inferred_tracts_list = []
total_true_tracts_list = []
# cutoff_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.9,0.95,0.96,0.97,0.98,0.985, 0.989, 0.99,0.995, 0.999, 0.9999, 0.99999, 0.999999]
cutoff_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.9,0.95, 0.99,0.995, 0.999]
# cutoff_list =  [0.999, 0.9999]
for cutoff in cutoff_list:
    total_inferred_tracts = 0
    total_true_tracts = 0
    true_positives = 0 
    for tree in tqdm.tqdm(range(100)):
#     for tree in tqdm.tqdm(range(100)):
        df = pd.read_csv(f'/scratch1/jcahoon/argml/data/distribution_{model_type}_validate/segments/{iter_num}/200kb_{iter_num}_{tree}.csv')
        for hap in range(56): # iterate the haplotypes
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
#             print( pred_sw.shape)
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
#             true_file = f'/scratch2/jcahoon/argml_data/{model}_50kb/bed_files/200kb_0_{tree}_{method}_{hap}_true.bed'
#             pred_file = f'/scratch2/jcahoon/argml_data/{model}_50kb/bed_files/200kb_0_{tree}_{method}_{hap}_pred.bed'
            true_file = f'/dev/shm/200kb_0_{tree}_{method}_{model_type}_{hap}_{iter_num}_true.bed'
            pred_file = f'/dev/shm/200kb_0_{tree}_{method}_{model_type}_{hap}_{iter_num}_pred.bed'
            true.sort(key=lambda x:x[0])
            pred.sort(key=lambda x:x[0])
#             print(true)
#             print(pred)
            with open(true_file, 'w') as o:
                for t in true:
                    o.write(f'1\t{t[0]}\t{t[1]}\n')
            with open(pred_file, 'w') as o:
                for t in pred:
                    o.write(f'1\t{t[0]}\t{t[1]}\n')
            it, tt, tp = cal_accuracy(true_file,pred_file)
            total_inferred_tracts += it
            total_true_tracts += tt
            true_positives += tp
    # the following is incorrect and should be reversed
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

with open(f'/project/lgarcia1_998/jcahoon/argml/pr_pickle/200kb_0_precision_recall_{model}_{method}_{model_type}_{iter_num}_{train_size*50}K.pickle', "wb") as f:
    pickle.dump([x,y, true_positives_list, total_inferred_tracts_list, total_true_tracts_list],f)