#!/usr/bin/env python2

import numpy as np
import pandas as pd
import random
import math
import sys
import vtktools

ntime = 989 # timesteps
uTot = np.array([])
for i in range(ntime):
    filename = '../data/small3DLSBU/LSBU_' + str(i) + '.vtu'
    print "Processing file " + filename 
    ug = vtktools.vtu(filename)
    ui = ug.GetScalarField('Tracer')
    # Create the 2D u matrix, [ntime,n]
    uTot = np.vstack([uTot, ui]) if uTot.size else ui 

print ug.GetFieldNames()
print uTot.shape

np.savetxt('../data/txt_data/background_state.txt', uTot)
# df.to_csv('../data/csv_data/background_state.csv')
