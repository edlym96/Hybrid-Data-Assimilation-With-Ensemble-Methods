#!/usr/bin/env python3
import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import time
import os

from numpy.linalg import inv
from numpy import linalg as LA

import math
from scipy.sparse.linalg import svds


import sys
#import vtktools
import argparse
#sys.path.append('fluidity-master')

from evaluate_DA_solution import evaluate_DA_solution

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

def build_DA_solution(xB_filepath, y_filepath, V_filepath, ntime = 989//2):


	xB = np.transpose(
		np.load(xB_filepath)['u'])  # Transpose from 989x100040 to 100040x989
	print("xB", xB.shape)

	n = xB.shape[0]
	y = np.transpose(
		np.load(y_filepath)['y'])  # Transpose from 989x100040 to 100040x989
	print("y", y.shape)

	# V = np.load('../data/matrix_prec_' + str(ntime) + '/matrixVprec' + str(trnc) + '.npz')['arr_0']
	V = np.load(V_filepath)['arr_0']
	lam = 0.1e-60
	R = lam * 0.9
	x0 = np.ones(n)

	Vin = np.linalg.pinv(V)
	print("Vin", Vin.shape)

	v0 = np.dot(Vin,x0)  # Take x0 from physical to reduced space by dotting with inverse of reduced V (This works because of the way the cost function is defined dx = Vdu)
	print("v0 ", v0.shape)
	VT = np.transpose(V)

	# Background state in obervation space. H is taken as identity matrix so a direct copy can be used.
	HxB = xB.copy()

	# Misfit calculated by subtracting
	d = np.subtract(y, HxB)
	#if localization:
		#build_Ch()
		#xB = np.multiply(Ch, xB)

	deltaxDA = np.zeros((n, ntime))

	t = time.time()
	for i in range(ntime):
		print("Processing timestep no: ", i)

		di = d[:,i].reshape([n,1])

		# Cost function J
		def J(v):
			# v = v.reshape((V.shape[1],ntime)) # need to reshape because optimize flattens input
			v = v.reshape((V.shape[1], 1))
			vT = np.transpose(v)
			vTv = np.dot(vT, v)
			Vv = np.dot(V, v).reshape([n,1])
			Jmis = np.subtract(Vv, di)
			invR = 1 / R
			JmisT = np.transpose(Jmis)
			RJmis = JmisT.copy()
			J1 = invR * np.dot(RJmis, Jmis)
			Jv = (vTv + J1) / 2
			# return LA.norm(Jv, 2)
			return Jv

		# Gradient of J
		# In this case, the adjoint operator, g is taken as identity
		def gradJ(v):
			# v = v.reshape((V.shape[1],ntime))
			v = v.reshape((V.shape[1], 1))
			Vv = np.dot(V, v).reshape([n,1])
			Jmis = np.subtract(Vv, di)
			invR = 1 / R
			# g1 = Jmis.copy()
			# g2 = np.dot(VT,g1)
			# gg2 = np.multiply(invR , g2)
			gg2 = np.multiply(invR, np.dot(VT, Jmis))  # VT[501x100040] Jmis[989x100040]
			ggJ = v + gg2
			# return ggJ.flatten()
			return ggJ


		res = minimize(J, v0, method='L-BFGS-B', jac=gradJ, options={'disp': False})
		vDA = np.reshape(res.x, (V.shape[1], 1))
		deltaxDA[:,i] = np.dot(V, vDA).flatten()  # take vDA from the reduced space back to x-space

	elapsed = time.time() - t
	print('elapsed', elapsed, 'seconds\n')
	xDA = xB + deltaxDA

	MSExDA = evaluate_DA_solution(xDA, xB, y)

	results_filename = os.path.basename(V_filepath)
	save_DA_solution(xDA, MSExDA, results_filename)

def save_DA_solution(xDA, MSE, filename):
	if not os.path.exists('../data/results'):
		os.makedirs('../data/results')
	print("Saving results to " + filename + "...")
	np.savez_compressed("../data/results/Results" + filename, xDA, result=MSE)
	

def build_Ch(pos_file, cutoff, rh):
	positions = pos_file.load(pos_file)['x']
	#Lh = max(positions) - min(positions)
	C = np.zeros([len(positions), len(positions)])
	for i in range(C.shape[0]):
		for j in range(i, C.shape[1]):
			s = abs(positions[i]-positions[j])
			if s <= cutoff/2:
				C[i,j] = C[j,i]= 1
			elif s >= cutoff:
				C[i,j] = C[j,i] = 0
			else:
				C[i,j] = C[j,i] = 0.5*(1+math.cos((2*math.pi*(s-cutoff/2))/cutoff))
	
	W,V = LA.eig(C)
	W = W[:,rh]
	V = np.sqrt(V[:rh,:rh])
	return np.matmul(W,V)

def arg_parser():
	parser = argparse.ArgumentParser(description='Calculate DA Solution')
	parser.add_argument('-xBp',
						'--xB_filepath',
						default="../data/converted_data/background_state.npz",
						help='provide file path for background state data'
						)
	parser.add_argument('-yp',
						'--y_filepath',
						default="../data/converted_data/observations.npz",
						help='provide file path for observation data'
						)
	parser.add_argument('-Vp',
						'--V_filepath',
						default="../data/matrix_prec_494/matrixVprec145.npz",
						help='provide file path for reduced error covariance matrix'
						)
	parser.add_argument('--ntime',
						default=494,
						help='number of timesteps'
						)
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = arg_parser()

	build_DA_solution(args.xB_filepath, args.y_filepath, args.V_filepath)


"""
ug.AddScalarField('uDA', xDA)
ug.Write('../VarDACode/Results/xDA-14Jun-3DTrac.vtu')

ug.AddScalarField('uM', xB)
ug.Write('../VarDACode/Results/xB-3DTrac.vtu')

ug.AddScalarField('v', y)
ug.Write('../VarDACode/Results/y-3DTrac.vtu')
"""