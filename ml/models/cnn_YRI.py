import torch
import torch.nn as nn
class CNN_YRI(nn.Module):
    def __init__(self):
        super().__init__()
        # YRI
        self.conv1_YRI = nn.Conv2d(1,4,5,2,2)
        self.pool_YRI = nn.MaxPool2d(2)
        self.conv2_YRI = nn.Conv2d(4,8,5,1,2)
        
        # CEU
        self.conv1_CEU = nn.Conv2d(1,4,5,2,2)
        self.pool_CEU = nn.MaxPool2d(2)
        self.conv2_CEU = nn.Conv2d(4,8,5,1,2)
        self.fc1 = nn.Linear(3136, 112) # for just CEU + YRI
    def _init_weights(self, module):
        if isinstance(module, nn.Linear):
            module.weight.data.normal_(mean=0.0, std=1.0)
            if module.bias is not None:
                module.bias.data.zero_()
    def forward(self, x):
        x[0] = self.pool_YRI(nn.functional.relu(self.conv1_YRI(x[0])))
        x[0] = self.pool_YRI(nn.functional.relu(self.conv2_YRI(x[0])))
        
        x[1] = self.pool_CEU(nn.functional.relu(self.conv1_CEU(x[1])))
        x[1] = self.pool_CEU(nn.functional.relu(self.conv2_CEU(x[1])))
        x = torch.cat((torch.flatten(x[0],1), torch.flatten(x[1],1)), 1)

        x = self.fc1(x)
        return x