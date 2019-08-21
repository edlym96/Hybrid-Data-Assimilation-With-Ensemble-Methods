#!/usr/bin/env python

import numpy as np
import os
import math
import sys
import vtktools
import argparse


def convert(filepath):
    ntime = 988 # timesteps (there are actually 989 but divide by 2 to split data evenly)
    features = 100040
    uTot = np.zeros([ntime//2, features])
    y = np.zeros([ntime//2, features])
    for i in range(ntime):
        filename = filepath + 'LSBU_' + str(i) + '.vtu'
        print("Processing file " + filename)
        ug = vtktools.vtu(filename)
        if i < ntime/2:
            ui = ug.GetScalarField('Tracer')
            # Create the 2D u matrix, [ntime,n]
            #uTot = np.vstack([uTot, ui]) if uTot.size else ui # model
            uTot[i%(ntime//2)] = ui
        else:
            yi = ug.GetScalarField('Tracer')
            #y = np.vstack([y,yi]) if y.size else yi # observations
            y[i%(ntime//2)] = yi

    pos = np.array(ug.GetLocations())

    print(ug.GetFieldNames())
    print(uTot.shape)
    print(y.shape)
    print(pos.shape)

    if not os.path.exists("../data"):
        os.mkdir("../data")
    if not os.path.exists("../data/converted_data"):
        os.mkdir("../data/converted_data")
    np.savez_compressed('../data/converted_data/background_state.npz', u=uTot)
    np.savez_compressed('../data/converted_data/observations.npz', y=y)
    np.savez_compressed('../data/converted_data/positions.npz', pos=pos)
    #np.savetxt('../data/txt_data/background_state.txt', uTot)
    # df.to_csv('../data/csv_data/background_state.csv')

def arg_parser():
    parser = argparse.ArgumentParser(description='Convert vtu to npz')
    parser.add_argument('-fp',
						'--filepath',
						default="../data/small3DLSBU/",
						help='provide file path data')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()

    convert(args.filepath)