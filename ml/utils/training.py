import torch
import numpy as np

def getMeanStd(test_list):
    mean = sum(test_list) / len(test_list) 
    variance = sum([((x - mean) ** 2) for x in test_list]) / len(test_list) 
    res = variance ** 0.5
    return mean, res

class W_BCEWithLogitsLoss(torch.nn.Module):
    def __init__(self, w_p = None, w_n = None):
        super(W_BCEWithLogitsLoss, self).__init__()

        self.w_p = w_p
        self.w_n = w_n

    def forward(self, logits, labels, epsilon = 1e-7):

        ps = torch.sigmoid(logits.squeeze()) 

        loss_pos = -1 * torch.mean(self.w_p * labels * torch.log(ps + epsilon))
        loss_neg = -1 * torch.mean(self.w_n * (1-labels) * torch.log((1-ps) + epsilon))

        loss = loss_pos + loss_neg

        return loss 
# dataseet loader
class eGRMDataset(torch.utils.data.Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor((X-Xmean)/Xstd)
        self.y = torch.tensor(y)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
#         return (self.X[None,:112,:112,idx],self.X[None,112:224,112:224,idx]), self.y[idx,:]
        return self.X[None,112:224,112:224,idx], self.y[idx,:]

# hamming distance loss
def Hamming_Loss(y_true, y_pred):
    temp=0
    for i in range(y_true.shape[0]):
        temp += np.size(y_true[i] == y_pred[i]) - np.count_nonzero(y_true[i] == y_pred[i])
    return temp/(y_true.shape[0] * y_true.shape[1])

def identity(x):
        return torch.tensor(x)