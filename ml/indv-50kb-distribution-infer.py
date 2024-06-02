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
os.chdir("~/argml/data/")
# EXAMPLE USAGE: python indv-50kb-distribution-growth-infer.py transformer growth 600000
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
# load files
if train_size == 1:
    url = f'~/arg-ml/data/distribution_{model_type}/tar_files/50kb_0_train.tar'
else:
    url = f'~/arg-ml/data/distribution_{model_type}/tar_files/50kb_{{0..{train_size-1}}}_train.tar'
dataset = (
    wds.WebDataset(url)
    .decode("l")
    .to_tuple("input.npy", "infer.input.npy", "output.npy")
    .map_tuple(preproc, identity)
)
trainloader = torch.utils.data.DataLoader(dataset, batch_size=128)
url = f'~/arg-ml/data/data/distribution_{model_type}/tar_files/50kb_{{0..3}}_validation.tar'
dataset = (
    wds.WebDataset(url)
    .decode("l")
    .to_tuple("input.npy", "infer.input.npy", "output.npy")
    .map_tuple(preproc, identity)
)
testloader = torch.utils.data.DataLoader(dataset, batch_size=128)

# initialize model
if model_name == "transformer":
    model = Transformer(d_model=112)
elif model_name == "cnn":
    model = CNN()
model.to(device)
trainable_params = sum(
	p.numel() for p in model.parameters() if p.requires_grad
)

print("Calculated weights: positive={}, negative={}".format(w_p, w_n))

# initialize loss function and optimizer
criterion = W_BCEWithLogitsLoss(w_p, w_n)
optimizer = torch.optim.Adam(model.parameters())

# training
with torch.enable_grad():
    n_epochs = 8
    for epoch in range(n_epochs):
        running_loss = 0.0
        y_pred = torch.empty(0)
        length = 0
        for i, (_, inputs, labels) in enumerate(trainloader,0):
            if model_name == "transformer":
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
        
# evaluation
y_pred = torch.empty(0)
y_out = torch.empty(0)
preds_y_list=[]
test_y_list =[]
with torch.no_grad():
    for i, (_, inputs, labels) in enumerate(testloader,0):
        if model_name == "transformer":
            inputs = inputs[:,:112,:112].to(device)
        else: 
            inputs = inputs[:,:,:112,:112].to(device)
        labels = labels[:, :112].to(device)
        outputs = model(inputs.float()) # prediction
        outputs = outputs.to("cpu")
        
        # saving outputs and labels
        pred_y = torch.round(torch.sigmoid(outputs))
        preds_y_list.append(outputs)
        y_pred=torch.cat((y_pred, pred_y))
        y_out=torch.cat((y_out, outputs))
        test_y_list.append(labels.to("cpu"))
y_test = np.concatenate(test_y_list, axis=0)
y_pred_t = y_pred.flatten()
y_test_t = y_test.flatten()

# DEBUG: print output
# print(torch.sigmoid(y_out))
# print(y_test)

# print test metrics
print("Test Metrics")
print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))

# save model
torch.save(model.state_dict(), f'trained_models/Distribution_{model_name}_{model_type}_50k_{train_size*50}K_infer.model')
