import argparse

import os
import torch
from torch import nn
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import numpy as np
import tqdm
import copy
import matplotlib.pyplot as plt
from torch.optim.lr_scheduler import ReduceLROnPlateau
import torchbnn as bnn

class EarlyStopping:
    """
    Early stopping handler for model training.

    This class implements early stopping functionality to stop model training
    when the validation loss stops improving.

    Args:
    - patience (int): Number of epochs to wait before stopping.
    - min_delta (float): Minimum change in loss to be considered an improvement.
    - verbose (bool): Whether to print early stopping messages.

    Methods:
    - __call__(val_loss): Check if early stopping criteria are met based on the validation loss.
    """

    def __init__(self, patience, min_delta, verbose=False):
        self.patience = patience
        self.min_delta = min_delta
        self.verbose = verbose
        self.counter = 0
        self.best_loss = None
        self.early_stop = False

    def __call__(self, val_loss):
        if self.best_loss is None:
            self.best_loss = val_loss
        elif val_loss > self.best_loss - self.min_delta:
            self.counter += 1
            if self.verbose:
                print(f"EarlyStopping counter: {self.counter} out of {self.patience}")
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_loss = val_loss
            self.counter = 0

def save_checkpoint(model, optimizer, epoch, loss, lr_scheduler, scaler, save_path):
    try:
        checkpoint = {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "lr_scheduler" : lr_scheduler,
            "loss": loss,
            "scaler_mean": scaler.mean_,
            "scaler_scale": scaler.scale_
        }
        torch.save(checkpoint, save_path)
    except IOError as e:
        print(f"Failed to save checkpoint to '{save_path}': {e}")

def load_checkpoint(model, optimizer, lr_scheduler, checkpoint_path, restart, scaler):
    if not os.path.exists(checkpoint_path):
        print("No checkpoint found. Starting from scratch.")
        return 0, False
    try:
        if restart:
            print("Restarting training, overwriting existing checkpoint.")
            return 0, False
        if map_location != None:
            checkpoint = torch.load(checkpoint_path, map_location=map_location)
        else:
            checkpoint = torch.load(checkpoint_path)
        model.load_state_dict(checkpoint["model_state_dict"])
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        lr_scheduler = checkpoint["lr_scheduler"]
        start_epoch = checkpoint["epoch"] + 1
        scaler.mean_ = checkpoint["scaler_mean"]
        scaler.scale_ = checkpoint["scaler_scale"]
        print(f"Checkpoint loaded. Resuming training from epoch {start_epoch}")
        return start_epoch, True
    except Exception as e:
        if "size mismatch" in str(e):
            error_msg = "Error: Checkpoint and current model do not match in size."
        else:
            error_msg = f"Error loading checkpoint from '{checkpoint_path}': {str(e)}"
        logging.error(error_msg)
        return 0, False

def get_options():
    description = 'Tests importance of classifiers from Pansim gridsearch'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python test_classifiers.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='Directory containing run files from gridsearch run')
    IO.add_argument('--outpref',
                    default = "output",
                    help='Output prefix. Default = "output"')
    IO.add_argument('--params',
                    required=True,
                    help='Comma separated list of predictors to estimate. Must match columns in infile')
    IO.add_argument('--hidden-size',
                    type=int,
                    default=512,
                    help='Size of hidden layers. Default = 512')
    IO.add_argument('--epochs',
                    type=int,
                    default=300,
                    help='Number of training epochs. Default = 100')
    IO.add_argument('--batch-size',
                    type=int,
                    default=10,
                    help='Batch size. Default = 10')
    parser.add_argument("--learning-rate", 
                        type=float, 
                        default=0.000001, 
                        help="Learning rate. Default = 0.000001")
    parser.add_argument("--lr-scheduler-factor",
                        type=float,
                        default=0.1, 
                        help="Factor by which the learning rate will be reduced by the learning rate scheduler")
    parser.add_argument("--lr-patience",
                        type=int, 
                        default=10, 
                        help="Patience for learning rate reduction")
    parser.add_argument("--weight-decay", 
                        type=float, 
                        default=1e-4, 
                        help="Weight decay for the optimizer")
    parser.add_argument("--early-stop-patience", 
                        type=int, 
                        default=10, 
                        help="Patience for early stopping")
    parser.add_argument("--min-delta", 
                        type=float, 
                        default=0.01, 
                        help="Minimum delta for early stopping")
    parser.add_argument("--bayesian", 
                        default=False,
                        action="store_true",
                        help="Fit Bayesian neural network.")
    parser.add_argument("--kl-weight", 
                        type=float, 
                        default=0.01, 
                        help="Weight for KL divergence in calculating loss.")
    parser.add_argument("--device", 
                        default=None, 
                        help="GPU device number if available. If not specified, will use all available Default = None")

    return parser.parse_args()


def get_device(device):
    if device is None:
        print("GPU not available, using cpu.")
        device = torch.device("cpu")
    else:
        if device != "cpu":
            device = torch.device("cuda:{}".format(device))
        else:
            device = torch.device("cpu")
    
    return device

class BayesianMultipleLinearRegression(torch.nn.Module):
    # Constructor
    def __init__(self, input_dim, output_dim, hidden_size):
        super(MultipleLinearRegression, self).__init__()
        self.linear = nn.Sequential(
            bnn.BayesLinear(prior_mu=0, prior_sigma=0.1, in_features=input_dim, out_features=hidden_size),
            nn.ReLU(),
            bnn.BayesLinear(prior_mu=0, prior_sigma=0.1, in_features=hidden_size, out_features=output_dim),
        )
        #self.linear = torch.nn.Linear(input_dim, output_dim)
    # Prediction
    def forward(self, x):
        y_pred = self.linear(x)
        return y_pred

class MultipleLinearRegression(torch.nn.Module):
    # Constructor
    def __init__(self, input_dim, output_dim, hidden_size):
        super(MultipleLinearRegression, self).__init__()
        self.linear = nn.Sequential(
            nn.Linear(input_dim, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, output_dim)
        )
        #self.linear = torch.nn.Linear(input_dim, output_dim)
    # Prediction
    def forward(self, x):
        y_pred = self.linear(x)
        return y_pred

def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref
    params = options.params.split(",")
    hidden_size = options.hidden_size
    n_epochs = options.epochs
    batch_size = options.batch_size
    learning_rate = options.learning_rate
    lr_scheduler_factor = options.lr_scheduler_factor
    weight_decay = options.weight_decay
    lr_patience = options.lr_patience
    early_stop_patience = options.early_stop_patience
    min_delta = options.min_delta
    kl_weight = options.kl_weight
    bayesian = options.bayesian

    # get device
    device = get_device(options.device)

    df = pd.read_csv(infile, header=0, sep="\t", )
    # set NAs to zero
    df = df.fillna(0)

    params = [col for col in params if col in df.columns]

    # Split data into X and Y parameters, convert bools to integers

    # generate training/test split
    train_df, test_df = train_test_split(df, train_size=0.7, shuffle=True)

    #print(f"train_df\n{train_df}")
    #print(f"test_df\n{test_df}")

    X_train_raw = train_df.drop(columns=params)*1
    X_test_raw = test_df.drop(columns=params)*1
    y_train = train_df[params]*1
    y_test = test_df[params]*1

    input_size = len(X_train_raw.columns)
    output_size = len(y_train.columns)

    X_train_raw = X_train_raw.to_numpy()

    # train scaler
    scaler = StandardScaler()
    scaler.fit(X_train_raw)
    
    X_train = torch.tensor(scaler.transform(X_train_raw), dtype=torch.float32)
    X_test = torch.tensor(scaler.transform(X_test_raw.to_numpy()), dtype=torch.float32)
    y_train = torch.tensor(y_train.to_numpy(), dtype=torch.float32)
    y_test = torch.tensor(y_test.to_numpy(), dtype=torch.float32)

    # print(f"X_train\n{X_train}\n{X_train.shape}")
    # print(f"X_test\n{X_test}\n{X_test.shape}")
    # print(f"y_train\n{y_train}\n{y_train.shape}")
    # print(f"y_test\n{y_test}\n{y_test.shape}")

    # initialise model
    if bayesian:
        model = BayesianMultipleLinearRegression(input_size, output_size, hidden_size)
    else:
        model = MultipleLinearRegression(input_size, output_size, hidden_size)

    # loss function and optimizer
    optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    # defining the loss criterion
    mse_loss = torch.nn.MSELoss()
    kl_loss = bnn.BKLLoss(reduction='mean', last_layer_only=False)
    #loss_fn = torch.nn.L1Loss()

    # lr scheduler
    lr_scheduler = ReduceLROnPlateau(optimizer, mode="min", factor=lr_scheduler_factor, patience=lr_patience) # taking big, then small steps

    # early stopping
    early_stopping = EarlyStopping(patience=early_stop_patience, min_delta=min_delta, verbose=True)

    # training parameters
    batch_start = torch.arange(0, len(X_train), batch_size)
    
    # Hold the best model
    best_loss = np.inf   # init to infinity
    best_weights = None
    history = []
    
    # training loop
    for epoch in range(n_epochs):
        model.train()
        with tqdm.tqdm(batch_start, unit="batch", mininterval=0, disable=False) as bar:
            bar.set_description(f"Epoch {epoch}")
            for start in bar:
                # take a batch
                X_batch = X_train[start:start+batch_size]
                y_batch = y_train[start:start+batch_size]
                # forward pass
                y_pred = model(X_batch)
                mse = mse_loss(y_pred, y_batch)
                loss = mse
                if bayesian:
                    kl = kl_loss(model)
                    loss += kl_weight*kl

                # print(f"loss\n{loss}")
                # if np.isnan(float(loss)) or np.isinf(float(loss)):
                #     print(f"y_pred\n{y_pred}")
                #     print(f"y_batch\n{y_batch}")
                #     fail
                # backward pass
                optimizer.zero_grad()
                loss.backward()
                # update weights
                optimizer.step()
                # print progress
                if bayesian:
                    bar.set_postfix({"mse": float(mse), "kl": float(kl.item()), "loss": float(loss.item())})
                else:
                    bar.set_postfix({"mse": float(mse), "loss": float(loss.item())})
        
        # evaluate accuracy at end of each epoch
        model.eval()
        y_pred = model(X_test)
        mse = mse_loss(y_pred, y_test)
        loss = mse
        if bayesian:
            kl = kl_loss(model)
            loss = mse + kl_weight*kl
            loss = float(loss.item())
            print(f"Epoch: {epoch} Test MSE: {float(mse)} KL: {float(kl.item())} Loss: {loss}")
        else:
            loss = float(loss.item())
            print(f"Epoch: {epoch} Test MSE: {float(mse)} Loss: {loss}")
        
        # archive loss
        history.append(loss)
        if loss < best_loss:
            best_loss = loss
            best_weights = copy.deepcopy(model.state_dict())
            save_checkpoint(model, optimizer, epoch, best_loss, lr_scheduler, scaler, outpref + ".chk")
        
        lr_scheduler.step(loss)
        early_stopping(loss)

        early_stop_tensor = torch.tensor(int(early_stopping.early_stop)).to(device)
        if early_stop_tensor.item() == 1:
            print("Early stopping triggered.", flush=True)
            break
            
    # restore model and return best accuracy
    #model.load_state_dict(best_weights)

    print("Loss: %.2f" % best_loss)
    #print("RMSE: %.2f" % np.sqrt(best_mse))
    plt.plot(history)
    plt.savefig(outpref + "_MSE.png")
    plt.close()

if __name__ == "__main__":
    main()