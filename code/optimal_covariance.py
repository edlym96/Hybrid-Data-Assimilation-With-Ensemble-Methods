#!/usr/bin/env python3
import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import time
import argparse
from scipy.linalg import inv
from numpy import linalg as LA
import time

import math
from scipy.sparse.linalg import svds

import os
import sys

# sys.path.append('fluidity-master')

# time steps
# ntime = 989
global ntime
""" Load the background state data from npz file """


def load_data(path):
    print("Loading data...")
    # xTot = np.loadtxt("../data/txt_data/background_state.txt")
    try:
        npzfile = np.load(path)
    except:
        print("Background state file not found. Please run convert_vtu.py")
        return
    xTot = npzfile['x']
    print("Shape of input data is ", xTot.shape)
    return xTot


"""BELOW IS FOR ENSEMBLE METHODS"""
"""get the background errors using ensemble"""


def ensemble_method(xTot, ensemble_size):
    t = time.time()
    print("Getting background error via ensemble...")
    ensemble = np.zeros([xTot.shape[1], ensemble_size])
    split_index = 0
    np.random.seed(7)
    # Split the data into groups to perform sampling
    for split in np.array_split(xTot, ensemble_size, axis=0):
        #  Get mean and std for sampling
        mean = np.mean(split, axis=0)
        std = np.std(split, axis=0)
        samples = np.zeros(mean.shape[0])
        for i in range(mean.shape[0]):
            samples[i] = np.random.normal(mean[i], std[i], 1)  # Sample for each feature
        ensemble[:, split_index] = samples
        split_index += 1
    ensemble_mean = np.expand_dims(np.mean(ensemble, axis=1), axis=1)  # Get ensemble mean
    Vens = (1 / math.sqrt(ensemble_size - 1)) * (ensemble - ensemble_mean)  # Get ensemble errors given by formula
    elapsed = time.time() - t
    print("Time taken for ensemble method " + str(elapsed))
    print("Ensemble is ", Vens.shape)
    return Vens


'''BELOW IS TSVD FORM OF SQUARE ROOT OF ERROR COVARIANCE'''


def tsvd_method(xTot, trnc):
    t = time.time()
    print("Getting background error via TSVD")
    err = np.absolute(xTot - np.mean(xTot, axis=0))
    print("err shape is ", err.shape)
    n = err.shape[1]

    V = np.transpose(err)  # V has the shape n x ntime

    Utrunc, strunc, Wtrunc = svds(V, k=trnc)
    Vtsvd = Utrunc.dot(np.diag(np.sqrt(strunc)))
    elapsed = time.time() - t
    print("Time taken for tsvd method is " + str(elapsed))
    print("V is ", V.shape)
    print("U is ", Utrunc.shape)
    print("s is ", strunc.shape)
    print("WT is ", Wtrunc.shape)
    print("U dot s is ", Vtsvd.shape)
    return Vtsvd


def save_covariance_matrices(V=None, Vens=None, ntime=494, trnc=145, ensemble_size=50, name=''):
    if not os.path.exists("../data/matrix_prec_" + str(ntime)):
        os.mkdir("../data/matrix_prec_" + str(ntime))
    print("Saving preconditioned background error covariances...")
    if V is not None:
        np.savez_compressed("../data/matrix_prec_" + str(ntime) + "/matrixVprec" + str(trnc) + name + ".npz", V=V)
    if Vens is not None:
        np.savez_compressed(
            "../data/matrix_prec_" + str(ntime) + "/matrixVensembleSplit" + str(ensemble_size) + name + ".npz", V=Vens)


def arg_parser():
    parser = argparse.ArgumentParser(description='Precondition Background Error Covariance Matrices')
    parser.add_argument('-fp',
                        '--filepath',
                        default="../data/converted_data/background_state.npz",
                        help='provide file path for background state data'
                        )
    parser.add_argument('-trnc',
                        default=145,
                        help='provide truncation parameter for tsvd'
                        )
    parser.add_argument('-ens_size',
                        default=40,
                        help='provide ensemble size'
                        )
    parser.add_argument('-ntime',
                        default=494,
                        help='number of time steps'
                        )
    parser.add_argument('-tsvd', action='store_true')
    parser.add_argument('-ens', action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()
    filepath = args.filepath
    trnc = int(args.trnc)
    ensemble_size = int(args.ens_size)

    ntime = args.ntime

    xTot = load_data(filepath)
    Vtsvd = None
    Vens = None

    if args.tsvd and not args.ens:
        Vtsvd = tsvd_method(xTot, trnc)
    elif args.ens and not args.tsvd:
        Vens = ensemble_method(xTot, ensemble_size)
    else:
        Vtsvd = tsvd_method(xTot, trnc)
        Vens = ensemble_method(xTot, ensemble_size)

    name = os.path.basename(filepath).replace("background_", '').replace(".npz", '')

    save_covariance_matrices(Vtsvd, Vens, ntime, trnc, ensemble_size, name)
