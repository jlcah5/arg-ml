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
from sklearn.metrics import accuracy_score, recall_score, f1_score, roc_auc_score, precision_score
from utils.training import  W_BCEWithLogitsLoss, Hamming_Loss, identity
from models.transformer import Transformer
from models.cnn import CNN
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
os.chdir("/scratch1/jcahoon/argml/data/detectSite/")

level = sys.argv[1]
model_name = sys.argv[2]
w_p = 0.9471964285714286
w_n = 0.05280357142857143
normalize = transforms.Normalize(
mean=[0.03383124259712062],
std=[0.21923093132895327])
preproc = transforms.Compose([
    transforms.ToTensor(),
    normalize
])


if level == "eval":
    url = f'/project/lgarcia1_998/jcahoon/SPS_50kb/50kb_{{1..19}}.tar'

    dataset = (
        wds.WebDataset(url)
        .decode("l")
    #     .rename(image="input", info="output")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    trainloader = torch.utils.data.DataLoader(dataset, batch_size=128)
    url = f'/project/lgarcia1_998/jcahoon/SPS_50kb/50kb_0.tar'
    dataset = (
        wds.WebDataset(url)
        .decode("l")
    #     .rename(image="input", info="output")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    testloader = torch.utils.data.DataLoader(dataset, batch_size=128)
elif level == "test":
    url = "/project/lgarcia1_998/jcahoon/SPS_50kb/50kb_0.tar"
    
    dataset = (
        wds.WebDataset(url)
        .decode("l")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    trainloader = torch.utils.data.DataLoader(dataset, batch_size=128)
    testloader = torch.utils.data.DataLoader(dataset, batch_size=128)
else:
    print("Invalid label. Choose test or eval.")
    sys.exit()
if model_name == "transformer":
    model = Transformer(d_model=112)
elif model_name == "cnn":
    model = CNN()
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)

print("Calculated weights: positive={}, negative={}".format(w_p, w_n))
criterion = W_BCEWithLogitsLoss(w_p, w_n)
optimizer = torch.optim.Adam(model.parameters())
# input_list=[]
# label_list=[]
with torch.enable_grad():
    n_epochs = 8
    if level == "test":
        n_epochs = 1
    for epoch in range(n_epochs):
        running_loss = 0.0
        y_pred = torch.empty(0)
        length = 0
        for i, (_, inputs, labels) in enumerate(trainloader,0):
            if model_name == "transformer":
#                 inputs_true = torch.squeeze(inputs[:,:,:112,:112]).to(device)
                inputs = inputs[:,:112,:112].to(device)
            else: 
                inputs = inputs[:,:,:112,:112].to(device)
            labels = labels[:, :112].to(device)
            optimizer.zero_grad()
            outputs = model(inputs.float())
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            outputs = outputs.to("cpu")
            pred_y = torch.round(torch.sigmoid(outputs))
            y_pred=torch.cat((y_pred, pred_y))
            length +=1
        print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss / length}')
    
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
print(torch.sigmoid(y_out))
print(y_test)
print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))
ofile = open(f'pr_curve/SPS_{model_name}_50k_{level}_infer_y_out.npy', "wb")
pickle.dump(torch.sigmoid(y_out), ofile)
ofile.close()

ofile = open(f'pr_curve/SPS_{model_name}_50k_{level}_infer_y_test.npy', "wb")
pickle.dump(y_test, ofile)
ofile.close()

torch.save(model.state_dict(), f'trained_models/SPS_{model_name}_50k_{level}_infer.model')

# print(y_pred)
# print(y_test)