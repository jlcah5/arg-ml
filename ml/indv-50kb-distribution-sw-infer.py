import torch
import torch.nn as nn
from torchvision import transforms
import webdataset as wds
import numpy as np
import pickle
import os
import sys
from sklearn.metrics import accuracy_score, recall_score, f1_score, precision_score
from utils.training import  W_BCEWithLogitsLoss, identity
from models.transformer import Transformer
from models.cnn import CNN

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
os.chdir("/project/lgarcia1_998/jcahoon/argml/ml")

# EXAMPLE USAGE: python indv-50kb-distribution-growth-sw-infer.py transformer growth 600000
# USER INPUT
model_name = sys.argv[1] # type of nn model: transformer OR cnn
model_type = sys.argv[2] # {growth, noGrowth}
train_size = int(int(sys.argv[3])/50000) # number of training samples, in increments of 50k
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
    
# itianlize model and load weights
if model_name == "transformer":
    model = Transformer(d_model=112)
elif model_name == "cnn":
    model = CNN()
model.load_state_dict(torch.load(f'trained_models/Distribution_{model_name}_{model_type}_50k_{level}_{train_size*50}K_infer.model'))
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)

# initialize lists to store metrics
acc_list = []
recall_list = []
prec_list = []

# iterate models and calculate performance for each evaluation set
# 4-13 are within distribution for 600k
# 100-109 are out of distribution for ALL models
for seed in [4,5,6,7,8,9,10,11,12,13, 100,101,102,103,104,105,106,107,108,109]:
    # load evluation set
    url = f'~/arg-ml/validation_data/{model_type}/200kb_{seed}.tar'
    dataset = (
        wds.WebDataset(url)
        .decode("l")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    testloader = torch.utils.data.DataLoader(dataset, batch_size=128)

    y_pred = torch.empty(0)
    y_out = torch.empty(0)
    preds_y_list=[]
    test_y_list =[]
    
    # prediction for evaluation set
    with torch.no_grad():
        for i, (_, inputs, labels) in enumerate(testloader,0):
            if model_name == "transformer":
                inputs = inputs[:,:112,:112].to(device)
            else: 
                inputs = inputs[:,:,:112,:112].to(device)
            labels = labels[:, :112].to(device)
            outputs = model(inputs.float())
            outputs = outputs.to("cpu")
            pred_y = torch.round(torch.sigmoid(outputs))
            preds_y_list.append(outputs)
            y_pred=torch.cat((y_pred, pred_y))
            y_out=torch.cat((y_out, outputs))
            test_y_list.append(labels.to("cpu"))
    y_test = np.concatenate(test_y_list, axis=0)
    y_pred_t = y_pred.flatten()
    y_test_t = y_test.flatten()
    print(f'Iter: {seed}----------------------------------------------')
    acc_list.append(accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
    recall_list.append( recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    prec_list.append(precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))
    
    # create files for pr_curve
    ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_out.npy', "wb")
    pickle.dump(torch.sigmoid(y_out), ofile)
    ofile.close()

    ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_test.npy', "wb")
    pickle.dump(y_test, ofile)
    ofile.close()
    
    # print metrics for set
    print(f'point based accuracy score:  {getMeanStd(acc_list)}')
    print(f'precision score:  {getMeanStd(prec_list)}')
    print(f'recall score:  {getMeanStd(recall_list)}')
# print average metrics
# in
print("----IN METRICS----")
print(f'point based accuracy score:  {getMeanStd(acc_list[:10])}')
print(f'precision score:  {getMeanStd(prec_list[:10])}')
print(f'recall score:  {getMeanStd(recall_list[:10])}')
# out
print("----OUT METRICS----")
print(f'point based accuracy score:  {getMeanStd(acc_list[10:])}')
print(f'precision score:  {getMeanStd(prec_list[10:])}')
print(f'recall score:  {getMeanStd(recall_list[10:])}')

