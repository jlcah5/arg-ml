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
from utils.training import  W_BCEWithLogitsLoss, Hamming_Loss, identity, getMeanStd
from models.transformer import Transformer
from models.cnn import CNN
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
os.chdir("/scratch1/jcahoon/argml/data/detectSite/")


model_name = sys.argv[1]
model_type = sys.argv[2]

w_p = 0.9471964285714286
w_n = 0.05280357142857143
normalize = transforms.Normalize(
mean=[0.033477144284946796],
std=[0.21640539250873986])
preproc = transforms.Compose([
    transforms.ToTensor(),
    normalize
])
if model_name == "transformer":
    model = Transformer(d_model=112)
elif model_name == "cnn":
    model = CNN()
if model_type == "true":
    model.load_state_dict(torch.load(f'trained_models/SPS_{model_name}_50k_eval.model'))
elif model_type == "transfer":
    model.load_state_dict(torch.load(f'trained_models/Gower_{model_name}_50k_eval_infer.model'))
    print("loading transfer")
else:
    model.load_state_dict(torch.load(f'trained_models/SPS_{model_name}_50k_eval_infer.model'))
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)
acc_list = []
recall_list = []
prec_list = []
for iter_num in range(10):

    url = f'/scratch2/jcahoon/argml_data/SPS_50kb/tarball/200kb_{iter_num}.tar'

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
    #                 inputs_true = torch.squeeze(inputs[:,:,:112,:112]).to(device)
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
    print(f'Iter: {iter_num}----------------------------------------------')
    acc_list.append(accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
    recall_list.append( recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    prec_list.append(precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
    print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))

#     print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
#     print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
#     print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
#     print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))
    ofile = open(f'pr_curve/SPS_{model_name}_50k_sw_{model_type}_{iter_num}_y_out.npy', "wb")
    pickle.dump(torch.sigmoid(y_out), ofile)
    ofile.close()

    ofile = open(f'pr_curve/SPS_{model_name}_50k_sw_{model_type}_{iter_num}_y_test.npy', "wb")
    pickle.dump(y_test, ofile)
    ofile.close()
print(f'point based accuracy score:  {getMeanStd(acc_list)}')
print(f'precision score:  {getMeanStd(prec_list)}')
print(f'recall score:  {getMeanStd(recall_list)}')

# with torch.no_grad():
#     for i, (inputs, _, labels) in enumerate(testloader,0):
#         if model_name == "transformer":
#                 inputs = torch.squeeze(inputs[:,:,:112,:112]).to(device)
# #                 inputs = inputs[:,:112,:112].to(device)
#         else: 
#             inputs = inputs[:,:,:112,:112].to(device)
#         labels = labels[:, :112].to(device)
#         outputs = model(inputs.float())
#         outputs = outputs.to("cpu")
#         pred_y = torch.round(torch.sigmoid(outputs))
#         preds_y_list.append(outputs)
#         y_pred=torch.cat((y_pred, pred_y))
#         y_out=torch.cat((y_out, outputs))
#         test_y_list.append(labels.to("cpu"))
# y_test = np.concatenate(test_y_list, axis=0)
# y_pred_t = y_pred.flatten()
# y_test_t = y_test.flatten()
# print(torch.sigmoid(y_out))
# print(y_test)
# print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
# print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
# print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
# print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))


# print(y_pred)
# print(y_test)