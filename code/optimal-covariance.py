#!/usr/bin/env python3
import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import time

from numpy.linalg import inv
from numpy import linalg as LA


import math
from scipy.sparse.linalg import svds

import os
import sys
#sys.path.append('fluidity-master')
#import vtktools

# time steps
#ntime = 989
ntime = 989//2
"""Get the background state data"""
"""
uvwTot = np.array([])
for i in range(ntime):
    filename = '../data/small3DLSBU/LSBU_'+str(i)+'.vtu'
    ug=vtktools.vtu(filename)
    ug.GetFieldNames()
    uvwVec=ug.GetScalarField('Tracer') # background state matrix
    n=len(uvwVec)
    uvwTot= np.append(uvwTot,uvwVec) # Get all background states for all times size of ntime x n

#  Get dimensions of uvwTot array
dimuvwTot = len(uvwTot)
"""

""" Load the background state data from npz file """
print("Loading data...")
#uTot = np.loadtxt("../data/txt_data/background_state.txt")
npzfile = np.load("../data/converted_data/background_state.npz")
uTot = npzfile['u']
print("Shape of input data is ", uTot.shape)

"""Get the background errors using mean"""
"""
m = np.array([])
m = np.zeros(dimuvwTot/ntime)
for j in range(n):
    for i in range(1,ntime+1):
        m[j] = np.copy(m[j] + uvwTot[j+(i-1)*n])
    m[j] = m[j]/ntime

err = np.array([])
err = np.zeros(dimuvwTot)

for j in range(n):
    for i in range(1,ntime+1):
        err[j+(i-1)*j] = abs(uvwTot[j+(i-1)*j]-m[j])
"""


"""BELOW IS FOR ENSEMBLE METHODS"""
"""get the background errors using ensemble"""
print("Getting background error via ensemble...")
ensemble_size = 50
mean = np.mean(uTot, axis = 0)
std = np.std(uTot, axis = 0)

np.random.seed(7)

ensemble = np.array([])
for i in range(mean.shape[0]):
#	print("Feature " + str(i))
	samples = np.random.normal(mean[i],std[i],ensemble_size)
	ensemble = np.vstack([ensemble, samples]) if ensemble.size else samples # ensemble is of size n x Nens	

ensemble_mean = np.expand_dims(np.mean(ensemble, axis = 1), axis = 1) # Get ensemble mean
Xens = (1/math.sqrt(ensemble_size-1))*(ensemble-ensemble_mean) # Get ensemble errors given by formula
print("Ensemble is ", Xens.shape)


'''BELOW IS TSVD FORM OF SQUARE ROOT OF ERROR COVARIANCE'''
print("Getting background error via TSVD")
err = np.absolute(uTot-np.mean(uTot, axis = 0))

print("err shape is ", err.shape)

n = err.shape[1]

#V =  np.transpose(np.reshape(err, (ntime,n)))
V = np.transpose(err) # V has the shape n x ntime (as stated in research paper)

'''
mid = int(np.floor(ntime/2))
U, s, W = svds(V, k=mid)
print('first eigh!')
st = np.sort(s)[::-1]
s1= st[0]
ref=math.sqrt(s1)
print "ref", ref
trnc=2
while st[trnc] >= ref and trnc < mid-1 :
        trnc=trnc+1
        print "trnc", trnc
        print "st[trnc]", st[trnc]
'''

trnc = 145


print('value of trnc')
print(trnc)


Utrunc, strunc, Wtrunc = svds(V, k=trnc)
X = Utrunc.dot(np.diag(np.sqrt(strunc)))

print("V is ", V.shape)
print("U is ", Utrunc.shape)
print("s is ", strunc.shape)
print("WT is ", Wtrunc.shape)
print("U dot s is ", X.shape)

if not os.path.exists("../data/matrix_prec_"+str(ntime)):
	os.mkdir("../data/matrix_prec_"+str(ntime))

print("Saving preconditioned background error covariances...")
np.savez_compressed("../data/matrix_prec_"+str(ntime)+"/matrixVprec"+str(trnc)+".npz", X)
np.savez_compressed("../data/matrix_prec_"+str(ntime)+"/matrixVensemble" + str(ensemble_size) + ".npz", Xens)


