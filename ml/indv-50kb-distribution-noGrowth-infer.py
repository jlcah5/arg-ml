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
from utils.training import  W_BCEWithLogitsLoss, Hamming_Loss, identity
from models.transformer import Transformer
from models.cnn import CNN
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)
os.chdir("/scratch1/jcahoon/argml/data/detectSite/")

level = sys.argv[1]
model_name = sys.argv[2]
model_type = "noGrowth"
train_size = int(int(sys.argv[3])/50000)
# train_type = sys.argv[3]
w_p = 1-0.06412648809523809
w_n = 0.06412648809523809
normalize = transforms.Normalize(
mean=[0.06262605265503221],
std=[0.20339860636469934])
preproc = transforms.Compose([
    transforms.ToTensor(),
    normalize
])


if level =="eval":
#     url = f'/scratch2/jcahoon/argml_data/Gower_50kb/tarball/50kb_{{1..19}}.tar'
    if train_size == 1:
        url = f'/scratch1/jcahoon/argml/data/distribution_{model_type}/tar_files/50kb_0_train.tar'
    else:
        url = f'/scratch1/jcahoon/argml/data/distribution_{model_type}/tar_files/50kb_{{0..{train_size-1}}}_train.tar'


    dataset = (
        wds.WebDataset(url)
#         .shuffle(1000)
        .decode("l")
    #     .rename(image="input", info="output")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    trainloader = torch.utils.data.DataLoader(dataset, batch_size=128)
#     url =  "/scratch2/jcahoon/argml_data/Gower_50kb/tarball/50kb_0.tar"
    url = f'/scratch1/jcahoon/argml/data/distribution_{model_type}/tar_files/50kb_{{0..3}}_validation.tar'
    dataset = (
        wds.WebDataset(url)
        .decode("l")
    #     .rename(image="input", info="output")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    testloader = torch.utils.data.DataLoader(dataset, batch_size=128)
    
elif level=="test":
    url = f'/scratch1/jcahoon/argml/data/distribution_{model_type}/tar_files/50kb_2_test.tar'
    
    dataset = (
        wds.WebDataset(url)
#         .shuffle(1000)
        .decode("l")
        .to_tuple("input.npy", "infer.input.npy", "output.npy")
        .map_tuple(preproc, identity)
    )
    trainloader = torch.utils.data.DataLoader(dataset, batch_size=128)
    testloader = torch.utils.data.DataLoader(dataset, batch_size=128)
#     evalloader = torch.utils.data.DataLoader(dataset, batch_size=128)

else:
    print("Invalid label. Choose test or train.")
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
#                 inputs = torch.squeeze(inputs[:,:,:112,:112]).to(device)
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
print(torch.sigmoid(y_out))
print(y_test)
print("Test Metrics")
print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))

# y_pred = torch.empty(0)
# y_out = torch.empty(0)
# preds_y_list=[]
# test_y_list =[]
# with torch.no_grad():
#     for i, (inputs, labels) in enumerate(evalloader,0):
#         if model_name == "transformer":
#                 inputs = torch.squeeze(inputs[:,:,:112,:112]).to(device)
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
# print("Eval Metrics")
# print("point based accuracy score: ", accuracy_score(y_pred=y_pred_t,y_true=y_test_t ))
# print("precision score: ", precision_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
# print("recall score: ", recall_score(y_pred=y_pred,y_true=y_test, average="samples", zero_division=0 ))
# print("f1 score: ", f1_score(y_pred=y_pred,y_true=y_test,  average="samples",  zero_division=0  ))
# ofile = open(f'pr_curve/Gower_{model_name}_50k_{level}_{train_type}_y_out.npy', "wb")
# pickle.dump(torch.sigmoid(y_out), ofile)
# ofile.close()

# ofile = open(f'pr_curve/Gower_{model_name}_50k_{level}_{train_type}_y_test.npy', "wb")
# pickle.dump(y_test, ofile)
# ofile.close()
ofile = open(f'pr_curve/Distribution_{model_name}_50k_{level}_{model_type}_{train_size*50}K_infer_y_out.npy', "wb")
pickle.dump(torch.sigmoid(y_out), ofile)
ofile.close()

ofile = open(f'pr_curve/Distribution_{model_name}_50k_{level}_{model_type}_{train_size*50}K_infer_y_test.npy', "wb")
pickle.dump(y_test, ofile)
ofile.close()
torch.save(model.state_dict(), f'trained_models/Distribution_{model_name}_{model_type}_50k_{level}_{train_size*50}K_infer.model')


# print(y_pred)
# print(y_test)