import math
import torch
import torch.nn as nn
import torch.nn.init as init 
import torch.nn.functional as F
import torch.optim as optim

def call_bn(bn, x):
    return bn(x)

class MLP(nn.Module):
    def __init__(self, input_size=2131, n_outputs=3, dropout_rate=0.5, top_bn=False):
        self.dropout_rate = dropout_rate
        self.top_bn = top_bn
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, 128)
        self.fc2 = nn.Linear(128, 32)
        self.fc3 = nn.Linear(32, n_outputs)

    def forward(self, x,):
        h=x
        h=self.fc1(h)
        h=F.leaky_relu(h, negative_slope=0.01)
        h = F.dropout2d(h, p=self.dropout_rate)

        h=self.fc2(h)
        h=F.leaky_relu(h, negative_slope=0.01)
        h = F.dropout2d(h, p=self.dropout_rate)

        logit=self.fc3(h)
        if self.top_bn:
            logit=call_bn(self.bn_c1, logit)
        return logit


