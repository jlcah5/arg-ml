"""
Description: make predictions for real data (note use python 3.7.7)
Usage: python3 indv-50kb-distribution-infer-real.py transformer growth PATH
Author: Jordan Cahoon, Sara Mathieson
Date: 6/21/24
"""

# python imports
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import seaborn as sns
import sys
import torch
from torch.utils.data import Dataset
from torchvision import transforms

# our imports
from models.cnn import CNN
from models.transformer import Transformer

# globals
NUM_HAPS = 112 # in target pop (224 in target + outgroup)
THRESH = 0.95  # threshold for high prob introgression

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

# USER INPUT
model_name = sys.argv[1] # type of nn model: transformer OR cnn
model_type = sys.argv[2] # {growth, noGrowth}
real_pop = sys.argv[3]   # path to real test data folder
out_prefix = sys.argv[4]   # path to output prefix (predictions as txt/pdf)

if model_type == "growth":
    # initialize parameters for loss function
    w_p = 1-0.05986458333333333 # proportion of positive haplotypes
    w_n = 0.05986458333333333 # proportion of negetive haplotypes

    # initilize normalization parameters or eGRM values
    normalize = transforms.Normalize(
    mean=[0.1087757664238845],
    std=[0.2459261757583878])
    preproc = transforms.Compose([
        transforms.ToTensor(),
        normalize
    ])
elif model_type == "noGrowth":
    # initialize parameters for loss function
    w_p = 1-0.06412648809523809 # proportion of positive haplotypes
    w_n = 0.06412648809523809 # proportion of negetive haplotypes
    
    # initilize normalization parameters or eGRM values
    normalize = transforms.Normalize(
    mean=[0.06262605265503221],
    std=[0.20339860636469934])
    preproc = transforms.Compose([
        transforms.ToTensor(),
        normalize
    ])
else:
    print("ERROR: model_type must be either growth or noGrowth")
    sys.exit()

class RealDataset(Dataset):

    def __init__(self, folder):
        unfiltered = os.listdir(folder)
        self.file_lst = []
        for f in unfiltered:
            if f[-3:] == "npy" and f[-6:-4] != "mu":
                self.file_lst.append(os.path.join(folder, f)) 

        # sort by length of filename, then by alpha/numeric
        self.file_lst.sort(key = lambda x: (len(x), x))

        self.file_lst = self.file_lst[:10] # try small subset

    def __len__(self):
        return len(self.file_lst)

    def __getitem__(self, idx):
        return np.load(self.file_lst[idx])

def parse_region(region_filename):
    # i.e. ends with GBR_chr1_2960000_3010000.npy
    tokens = region_filename.split("/")[-1].split("_")
    chrom = tokens[1][3:]
    start = tokens[2]
    end = tokens[3].split(".")[0]
    print(chrom, start, end)
    return chrom, start, end

dataset = RealDataset(real_pop)
testloader = torch.utils.data.DataLoader(dataset, batch_size=128)

# initialize model and load weights
if model_name == "transformer":
    model = Transformer(d_model=NUM_HAPS)
elif model_name == "cnn":
    model = CNN()
model.load_state_dict(torch.load(f'trained_models/Distribution_{model_name}_{model_type}_50k_eval_600K_infer.model'))
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)

y_pred = torch.empty(0) # rounded
y_out = torch.empty(0)  # probability
preds_y_list=[]
#test_y_list =[]

# prediction for real dataset
with torch.no_grad():
    for i, inputs in enumerate(testloader):
        if model_name == "transformer":
            inputs = inputs[:,:NUM_HAPS,:NUM_HAPS].to(device)
        else: 
            inputs = inputs[:,:,:NUM_HAPS,:NUM_HAPS].to(device)
        print(inputs.shape)
        #labels = labels[:, :NUM_HAPS].to(device)
        outputs = model(inputs.float())
        outputs = outputs.to("cpu")
        pred_y = torch.round(torch.sigmoid(outputs))
        preds_y_list.append(outputs)
        y_pred=torch.cat((y_pred, pred_y))
        y_out=torch.cat((y_out, torch.sigmoid(outputs)))
        #test_y_list.append(labels.to("cpu"))
#y_test = np.concatenate(test_y_list, axis=0)
y_pred_t = y_pred.flatten().numpy()
y_out_t = y_out.flatten().numpy()
num_pred = y_out_t.shape[0]
num_high = (y_out_t >= THRESH).sum()
print("frac high", num_high, "/", num_pred, num_high/num_pred)
sns.displot(y_out_t)#, x="prob introgression")
plt.xlabel("prob introgression")
plt.title(model_name + ", " + model_type + ", " + real_pop[-4:-1])
plt.savefig(out_prefix + ".pdf", bbox_inches="tight")

# save to a file
pred_file = open(out_prefix + ".txt", 'w')
file_lst = dataset.file_lst
assert len(file_lst)*NUM_HAPS == len(y_out_t)

i = 0
for f in file_lst:
    chrom, start, end = parse_region(f)
    for hap in range(NUM_HAPS):
        pred = y_out_t[i]
        line = "\t".join([chrom, start, end, str(pred), str(hap)])
        pred_file.write(line + "\n")
        i += 1

pred_file.close()

#y_test_t = y_test.flatten()
#print(f'Iter: {seed}----------------------------------------------')
#acc_list.append(accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
#recall_list.append( recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
#prec_list.append(precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
#print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))

# create files for pr_curve
#ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_out.npy', "wb")
#pickle.dump(torch.sigmoid(y_out), ofile)
#ofile.close()

#ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_test.npy', "wb")
#pickle.dump(y_test, ofile)
#ofile.close()
