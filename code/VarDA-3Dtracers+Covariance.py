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


import sys
sys.path.append('fluidity-master')
#import vtktools


#ntime = 989 # Set the number of time steps
ntime = 989//2 # Set the number of time steps
trnc = 145 #501 # Set the truncation parameter

m = trnc

#  V is the truncated form of square root of the background error covariance, B
#V = np.loadtxt('/.../matrixV'+str(m)+'-velocity.txt', usecols=range(m))
#V = np.load('../data/matrix_prec_' + str(ntime) + '/matrixVprec' + str(trnc) + '.npz')['arr_0']
V = np.load('../data/matrix_prec_' + str(ntime) + '/matrixVensemble60.npz')['arr_0']

print("V", V.shape)

"""
# Read .vtu file
ugg=vtktools.vtu('/homes/rarcucci/4DVAR-ROX/VarDACode/small3DCase/LSBU_100.vtu')
# Get the field names for the vtu file
ugg.GetFieldNames()
# Get tracer observations
uvwVecobs = ugg.GetScalarField('Tracer')

ug=vtktools.vtu('/homes/rarcucci/4DVAR-ROX/VarDACode/small3DCase/LSBU_10.vtu')
ug.GetFieldNames()
# Get tracer values
uvwVec = ug.GetScalarField('Tracer')

# Get coordinates
pos=ug.GetLocations()
# Get z-axis values
z=pos[:,2]
"""

lam = 0.1e-60

#n = len(uvwVec)

# Background state matrix
#xB = np.expand_dims(np.load('../data/converted_data/background_state.npz')['u'][10], axis = 0)
xB = np.transpose(np.load('../data/converted_data/background_state.npz')['u']) # Transpose from 989x100040 to 100040x989
n = xB.shape[0]

print("xB", xB.shape)

# Observations
#y = np.expand_dims(np.load('../data/converted_data/observations.npz')['y'][10], axis = 0)
y = np.transpose(np.load('../data/converted_data/observations.npz')['y'])  # Transpose from 989x100040 to 100040x989 
print("y", y.shape)

# Observation error covariance matrix
R = lam * 0.9

# Initial conditions set to a vector of ones
#x0 = np.ones(n)
x0 = np.ones([n,ntime])

# Get inverse of V
Vin = np.linalg.pinv(V)

print("Vin", Vin.shape)

"""v0 is u0 in research paper. """
v0 = np.dot(Vin,x0)  # Take x0 from physical to reduced space by dotting with inverse of reduced V (This works because of the way the cost function is defined dx = Vdu)
# Get transpose of V
VT = np.transpose(V)

# Background state in obervation space. H is taken as identity matrix so a direct copy can be used.
HxB = xB.copy()

# Misfit calculated by subtracting
d = np.subtract(y,HxB) # Transpose from 989x100040 to 100040x989

print("v0 ", v0.shape)
print("d", d.shape)

"""Right now, cost function uses approximate form of Background error covariance matrix, B
	Need to change it from using B=VVT to Pf where Pf is an ensemble of error covariances.
	Should build ensemble Pf using u from optimal-covariance.py file.
	Find C to perform localisation
	Change Cost function and grad
"""

# Need to structure the cost function to take in 1d input and reshape it back into matrix

# Cost function J
def J(v):
	v = v.reshape((V.shape[1],ntime)) # need to reshape because optimize flattens input
	#v = v.reshape((V.shape[1],1)) 
	vT = np.transpose(v)
	vTv = np.dot(vT,v)
	Vv = np.dot(V,v)
	Jmis = np.subtract(Vv,d)
	invR = 1/R
	JmisT = np.transpose(Jmis)
	RJmis = JmisT.copy()
	J1 = invR*np.dot(RJmis,Jmis)
	Jv = (vTv + J1) / 2
	return LA.norm(Jv,2)


# Gradient of J
# In this case, the adjoint operator, g is taken as identity
def gradJ(v):
	v = v.reshape((V.shape[1],ntime))
	#v = v.reshape((V.shape[1],1)) 
	Vv = np.dot(V,v)
	Jmis = np.subtract(Vv,d)
	invR = 1/R
	#g1 = Jmis.copy()
	VT = np.transpose(V)
	#g2 = np.dot(VT,g1)
	#gg2 = np.multiply(invR , g2)
	gg2 = np.multiply(invR,np.dot(VT,Jmis)) #VT[501x100040] Jmis[989x100040]
	ggJ = v + gg2
	return ggJ.flatten()

# Compute the minimum
t = time.time()

res = minimize(J, v0, method='L-BFGS-B', jac=gradJ, options={'disp': True})

vDA = np.array([])
vDA = np.reshape(res.x,(V.shape[1],ntime))

deltaxDA = np.dot(V,vDA) # take vDA from the reduced space back to x-space
xDA = xB + deltaxDA

elapsed = time.time() - t
print('elapsed', elapsed, '\n')

errxB = y - xB
MSExb = LA.norm(errxB, 2)/LA.norm(y, 2)
print('L2 norm of the background error' , MSExb , '\n')

errxDA = y - xDA
MSExDA = LA.norm(errxDA, 2)/LA.norm(y, 2)
print('L2 norm of the error in DA solution' , MSExDA , '\n')

"""
ug.AddScalarField('uDA', xDA)
ug.Write('../VarDACode/Results/xDA-14Jun-3DTrac.vtu')

ug.AddScalarField('uM', xB)
ug.Write('../VarDACode/Results/xB-3DTrac.vtu')

ug.AddScalarField('v', y)
ug.Write('../VarDACode/Results/y-3DTrac.vtu')
"""
