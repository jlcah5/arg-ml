import torch
import torch.nn as nn
import math

class PositionalEncoding(nn.Module):
    """
    https://pytorch.org/tutorials/beginner/transformer_tutorial.html
    """

    def __init__(self, d_model, vocab_size=5000, dropout=0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        pe = torch.zeros(vocab_size, d_model)
        position = torch.arange(0, vocab_size, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float()
            * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer("pe", pe)

    def forward(self, x):
#         print(x.shape)
        x = x + self.pe[:, : x.size(1), :]
        return self.dropout(x)
    
# transformer model
class Transformer(nn.Module):
    """
    Text classifier based on a pytorch TransformerEncoder: https://n8henrie.com/2021/08/writing-a-transformer-classifier-in-pytorch/
    """

    def __init__(
        self,
        d_model,
        nhead=8,
        dim_feedforward=8, # increase this ??
        num_layers=1,
        dropout=0.1,
        activation="relu",
        classifier_dropout=0.1,
        encoding=False
    ):

        super().__init__()
        
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first = True
        )
        self.pos_encoder = PositionalEncoding(
            d_model=d_model,
            dropout=dropout,
            vocab_size=128,
        )
        self.transformer_encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_layers
        )
        self.classifier = nn.Linear(d_model, 112) # could change this to a different classifer ??
        self.d_model = d_model
        self.encoding = encoding

    def forward(self, x):
#         x = self.emb(x) * math.sqrt(self.d_model)
        if self.encoding == True:
            x = self.pos_encoder(x)
        x = self.transformer_encoder(x)
        x = x.mean(dim=1) # could try max 
        x = self.classifier(x)
        return x