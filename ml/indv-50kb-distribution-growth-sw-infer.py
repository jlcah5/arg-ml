import torch
import torch.nn as nn
from torch.utils.data import IterableDataset
from torchvision import transforms
import webdataset as wds
import numpy as np
import pickle
import os
import sys
import math
from sklearn.model_selection import train_test_split
from torch.utils.data.dataloader import default_collate
from sklearn.metrics import accuracy_score, recall_score, f1_score, roc_auc_score, precision_score
from utils.training import  W_BCEWithLogitsLoss, Hamming_Loss, identity,  getMeanStd
from models.transformer import Transformer
from models.cnn import CNN
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
# os.chdir("/scratch1/jcahoon/argml/data/detectSite/")
os.chdir("/project/lgarcia1_998/jcahoon/argml/ml")

level = "eval"
model_name = sys.argv[1]
model_type = "growth"
train_size = int(int(sys.argv[2])/50000)
# train_type = sys.argv[3]
w_p = 1-0.05986458333333333
w_n = 0.05986458333333333
normalize = transforms.Normalize(
mean=[0.1087757664238845],
std=[0.2459261757583878])
preproc = transforms.Compose([
    transforms.ToTensor(),
    normalize
])

if model_name == "transformer":
    model = Transformer(d_model=112)
elif model_name == "cnn":
    model = CNN()
model.load_state_dict(torch.load(f'trained_models/Distribution_{model_name}_{model_type}_50k_{level}_{train_size*50}K_infer.model'))
# if model_type == "true":
#     model.load_state_dict(torch.load(f'trained_models/Gower_{model_name}_50k_eval.model'))
# elif model_type == "transfer":
#     model.load_state_dict(torch.load(f'trained_models/SPS_{model_name}_50k_eval_infer.model'))
#     print("loading transfer")
# else:
#     model.load_state_dict(torch.load(f'trained_models/Gower_{model_name}_50k_eval_infer.model'))
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)
acc_list = []
recall_list = []
prec_list = []
# for seed in [80,81,82,83,84,85,86,87,88,89,100,101,102,103,104,105,106,107,108,109]: 
# for seed in [100,101,102,103,104,105,106,107,108,109]: 
# for seed in [4,5,6,7,8,9,10,11,12,13]: 
for seed in [11,12]: 
    url = f'/scratch1/jcahoon/argml/data/distribution_growth_validate/tar_filelist/200kb_{seed}.tar'

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
    with torch.no_grad():
        for i, (_, inputs, labels) in enumerate(testloader,0):
            if model_name == "transformer":
    #                 inputs = torch.squeeze(inputs[:,:,:112,:112]).to(device)
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
    #     print(torch.sigmoid(y_out))
    #     print(y_test)
    print(f'Iter: {seed}----------------------------------------------')
    acc_list.append(accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
    recall_list.append( recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    prec_list.append(precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))

    #     print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
    #     print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    #     print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    #     print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))
    ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_out.npy', "wb")
    pickle.dump(torch.sigmoid(y_out), ofile)
    ofile.close()

    ofile = open(f'pr_curve/Distribution_{model_name}_50k_sw_{model_type}_{seed}_{train_size*50}K_y_test.npy', "wb")
    pickle.dump(y_test, ofile)
    ofile.close()

    print(f'point based accuracy score:  {getMeanStd(acc_list)}')
    print(f'precision score:  {getMeanStd(prec_list)}')
    print(f'recall score:  {getMeanStd(recall_list)}')
# in
print("----IN METRICS----")
print(f'point based accuracy score:  {getMeanStd(acc_list[:10])}')
print(f'precision score:  {getMeanStd(prec_list[:10])}')
print(f'recall score:  {getMeanStd(recall_list[:10])}')
# out
# print("----OUT METRICS----")
# print(f'point based accuracy score:  {getMeanStd(acc_list[10:])}')
# print(f'precision score:  {getMeanStd(prec_list[10:])}')
# print(f'recall score:  {getMeanStd(recall_list[10:])}')

# print(y_pred)
# print(y_test)