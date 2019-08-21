#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize
import random
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os

from numpy.linalg import inv
from scipy import linalg as LA

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

def build_DA_solution(xB_filepath, y_filepath, V_filepath, pos_filepath, ntime = 989//2, localisation = True):


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
	# x0 = np.ones(n)
	x0 = np.zeros(n)

	# Background state in obervation space. H is taken as identity matrix so a direct copy can be used.
	HxB = xB.copy()

	# Misfit calculated by subtracting
	d = np.subtract(y, HxB)

	if localisation:
		x_pos = np.load(pos_filepath)['pos'][:,0]
		y_pos = np.load(pos_filepath)['pos'][:,1]
		z_pos = np.load(pos_filepath)['pos'][:,2]
		Ch = localise_h(x_pos, y_pos, 250) #Ch is (100000 , 50)
		Cv = localise_v(z_pos, 10)
		V_new = np.zeros([V.shape[0], V.shape[1]*Ch.shape[1]*Cv.shape[1]])
		for i in range(V.shape[1]):
			tmp = np.tile(V[:,i],(Ch.shape[1],1)).transpose() #Shape (100000,50)
			tmp = np.multiply(Ch, tmp)
			for j in range(tmp.shape[1]):
				tmp2 = np.tile(tmp[:,j],(Cv.shape[1],1)).transpose() #Shape (100000,40)
				tmp2 = np.multiply(Cv,tmp2)
				V_new[:,(i*Ch.shape[1]+j)*Cv.shape[1]:(i*Ch.shape[1]+j+1)*Cv.shape[1]] = tmp2

		V = V_new

	Vin = np.linalg.pinv(V)

	print("Vin", Vin.shape)

	v0 = np.dot(Vin,
				x0)  # Take x0 from physical to reduced space by dotting with inverse of reduced V (This works because of the way the cost function is defined dx = Vdu)
	print("v0 ", v0.shape)
	VT = np.transpose(V)

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
	save_DA_solution(xDA, MSExDA, results_filename, localisation)


def save_DA_solution(xDA, MSE, filename, localisation):
	if not os.path.exists('../data/results'):
		os.makedirs('../data/results')
	print("Saving results to " + filename + "...")
	path = "../data/results/Results"
	path += "Localisation" if localisation else ''
	path += filename
	np.savez_compressed(path, xDA=xDA, result=MSE)


def localise_h(x_positions, y_positions, cutoff, rh=10):
	#Lh = max(positions) - min(positions)
	# C = np.zeros([len(positions), len(positions)])
	# TODO: MOVE CONSTRUCTION OF LOCALISATION MATRIX TO OPTIMAL COVARIANCE
	# downsampled_x = x_positions[::scale]
	# downsampled_y = y_positions[::scale]

	# positions = np.stack((downsampled_x,downsampled_y), axis=1)
	positions = np.stack((x_positions, y_positions), axis=1)
	C = np.zeros((positions.shape[0],positions.shape[0]))
	# t = time.time()
	for i in range(positions.shape[0]):
		coord = positions[i]
		for j in range(i,positions.shape[0]):
			s = LA.norm(positions[j]-coord, ord=2)
			if s <= cutoff / 2:
				C[i,j] = C[j,i]= 1
			elif s > cutoff:
				C[i,j] = C[j,i] = 0
			else:
				C[i,j] = C[j,i] = 0.5*(1+math.cos((2*math.pi*(s-cutoff/2))/cutoff))


	W, V = LA.eigh(C, eigvals=(C.shape[0]-rh, C.shape[0]-1))
	# print(C)
	# print(W)
	# print(V)
	# print(W.shape, V.shape)
	# print(np.matmul(V, np.diag(W)).shape)
	idx = W.argsort()[::-1]
	# W = np.sqrt(W[idx])
	W = W[idx]
	V = V[:, idx]
	# W = np.sqrt(W)
	# print(time.time()-t)
	# # print(np.matmul(W, V).shape)
	return np.matmul(V, np.diag(W))


def localise_v(z_positions, scale, rv=10):
	C = np.zeros((z_positions.shape[0], z_positions.shape[0]))
	for i in range(z_positions.shape[0]):
		coord = z_positions[i]
		for j in range(i,z_positions.shape[0]):
			dz = abs(z_positions[j]-coord)
			C[i, j] = C[j, i] = 1/(1+(dz/scale)**2)
	W, V = LA.eigh(C, eigvals=(C.shape[0] - rv, C.shape[0] - 1))
	idx = W.argsort()[::-1]
	# W = np.sqrt(W[idx])
	W = W[idx]
	V = V[:, idx]
	return np.matmul(V,np.diag(W))

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
	parser.add_argument('-posp',
						'--pos_filepath',
						default="../data/converted_data/positions.npz",
						help='provide file path for position matrix'
						)
	parser.add_argument('--ntime',
						default=494,
						help='number of timesteps'
						)
	parser.add_argument('--localisation',
						'-local',
						action='store_false',
						help='Turn on localisation'
						)
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = arg_parser()

	build_DA_solution(args.xB_filepath, args.y_filepath, args.V_filepath, args.pos_filepath, args.ntime, args.localisation)


"""
ug.AddScalarField('uDA', xDA)
ug.Write('../VarDACode/Results/xDA-14Jun-3DTrac.vtu')

ug.AddScalarField('uM', xB)
ug.Write('../VarDACode/Results/xB-3DTrac.vtu')

ug.AddScalarField('v', y)
ug.Write('../VarDACode/Results/y-3DTrac.vtu')
"""
