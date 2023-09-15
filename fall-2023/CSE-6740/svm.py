# implement support vector classifier
# Modified from https://github.com/kazuto1011/svm-pytorch/blob/master/main.py (Kazuto Nakashima) MIT License


import argparse
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.datasets import make_blobs
from torch.utils.data import Dataset, DataLoader

def create_model():
    return nn.Linear(2,1).to(args.device)

def create_data():
    X, Y = make_blobs(n_samples=100, n_features=2, centers=2, random_state=2)
    Y[Y == 0] = -1
    return torch.tensor(X, dtype=torch.float32), torch.tensor(Y, dtype=torch.float32)

def create_dataset(X, Y):
    class dataset(Dataset):
        def __init__(self, X, Y):
            self.X = X
            self.Y = Y
        def __len__(self):
            return len(self.Y)
        def __getitem__(self, idx):
            return self.X[idx], self.Y[idx]
    return dataset(X,Y)
    

def train(X, Y, model, args):
    optimizer = optim.SGD(model.parameters(), lr=args.lr)
    trainset = create_dataset(X, Y)
    trainloader = DataLoader(trainset, batch_size=args.batchsize, shuffle=True)
    N = len(trainset)
    model.train()
    for epoch in range(args.epoch):
        sum_loss = 0.0
        for t, (x, y) in enumerate(trainloader):
            x_g = x.to(args.device)
            y_g = y.to(args.device)
            output = model(x_g).squeeze()
            weight = model.weight.squeeze()

            loss = torch.mean(torch.clamp(1 - y_g * output, min=0))
            loss += args.c * (weight.t() @ weight) / 2.0

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            sum_loss += float(loss)

        print("Epoch: {:4d}\tloss: {}".format(epoch, sum_loss / N))


def visualize(X, Y, model):
    W = model.weight.squeeze().detach().cpu().numpy()
    b = model.bias.squeeze().detach().cpu().numpy()

    delta = 0.001
    x1 = np.linspace(X[:, 0].min(), X[:, 0].max(), 200)
    x2 = -b / W[1] - W[0] / W[1] * x1
    x2_p = (1-b) / W[1] - W[0] / W[1] * x1 
    x2_m = (-1-b) / W[1] - W[0] / W[1] * x1
    
    x1_data = X[:, 0].detach().cpu().numpy()
    x2_data = X[:, 1].detach().cpu().numpy()
    y_data = Y.detach().cpu().numpy()

    x1_c1 = x1_data[np.where(y_data == 1)]
    x1_c2 = x1_data[np.where(y_data == -1)]
    x2_c1 = x2_data[np.where(y_data == 1)]
    x2_c2 = x2_data[np.where(y_data == -1)]

    plt.figure(figsize=(10, 10))
    plt.xlim([x1_data.min() + delta, x1_data.max() - delta])
    plt.ylim([x2_data.min() + delta, x2_data.max() - delta])
    plt.plot(x1_c1, x2_c1, "P", c="b", ms=10)
    plt.plot(x1_c2, x2_c2, "o", c="r", ms=10)
    plt.plot(x1, x2, "k", lw=3)
    plt.plot(x1, x2_p, "b--", lw=3)
    plt.plot(x1, x2_m, "r--", lw=3)
    plt.tight_layout()
    plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=0.01)
    parser.add_argument("--lr", type=float, default=0.1)
    parser.add_argument("--batchsize", type=int, default=5)
    parser.add_argument("--epoch", type=int, default=10)
    parser.add_argument("--device", default="mps", choices=["cpu", "cuda", "mps"])
    args = parser.parse_args()
    if (args.device == "mps" and not torch.backends.mps.is_available()) or (args.device == "cuda" and not torch.cuda.is_avaliable()):
        args.device = "cpu"
    args.device = torch.device(args.device)

    print(args)

    X, Y = create_data()
    model = create_model()
    train(X, Y, model, args)
    visualize(X, Y, model)
