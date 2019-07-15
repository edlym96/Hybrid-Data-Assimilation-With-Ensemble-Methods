#!/usr/bin/env python2

import numpy as np
import pandas as pd
import random
import math
import sys
import vtktools

ntime = 988 # timesteps (there are actually 989 but divide by 2 to split data evenly)
uTot = np.array([])
y = np.array([])
for i in range(ntime):
    filename = '../data/small3DLSBU/LSBU_' + str(i) + '.vtu'
    print "Processing file " + filename 
    ug = vtktools.vtu(filename)
    if i < ntime/2:
    	ui = ug.GetScalarField('Tracer')
    	# Create the 2D u matrix, [ntime,n]
    	uTot = np.vstack([uTot, ui]) if uTot.size else ui # model
    else:
    	yi = ug.GetScalarField('Tracer')
    	y = np.vstack([y,yi]) if y.size else yi # observations

print ug.GetFieldNames()
print uTot.shape
print y.shape

np.savez_compressed('../data/converted_data/background_state.npz', u=uTot)
np.savez_compressed('../data/converted_data/observations.npz', y=y)
#np.savetxt('../data/txt_data/background_state.txt', uTot)
# df.to_csv('../data/csv_data/background_state.csv')
