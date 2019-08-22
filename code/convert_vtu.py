#!/usr/bin/env python

import numpy as np
import os
import math
import sys
import vtktools
import argparse
from scipy import linalg as LA

def convert(filepath, field_name='Tracer'):
    ntime = 988 # timesteps (there are actually 989 but divide by 2 to split data evenly)
    features = 100040
    uTot = np.zeros([ntime//2, features])
    y = np.zeros([ntime//2, features])
    for i in range(ntime):
        filename = filepath + 'LSBU_' + str(i) + '.vtu'
        print("Processing file " + filename)
        ug = vtktools.vtu(filename)
        if i < ntime/2:
            if field_name != 'Velocity':
                ui = ug.GetScalarField(field_name)

                # Create the 2D u matrix, [ntime,n]
                #uTot = np.vstack([uTot, ui]) if uTot.size else ui # model
                uTot[i%(ntime//2)] = ui
            else:
                ui = np.array(ug.GetVectorField(field_name))
                ui = LA.norm(ui, 2, axis=1)

                uTot[i % (ntime // 2)] = ui
        else:
            if field_name != 'Velocity':
                yi = ug.GetScalarField(field_name)
                #y = np.vstack([y,yi]) if y.size else yi # observations
                y[i%(ntime//2)] = yi
            else:
                yi = ug.GetVectorField(field_name)
                yi = LA.norm(yi, 2, axis=1)

                y[i % (ntime // 2)] = yi


    print(ug.GetFieldNames())
    print(uTot.shape)
    print(y.shape)

    if not os.path.exists("../data"):
        os.mkdir("../data")
    if not os.path.exists("../data/converted_data"):
        os.mkdir("../data/converted_data")
    if field_name == 'Tracer':
        back_path = '../data/converted_data/background_state.npz'
        obs_path = '../data/converted_data/observations.npz'
    elif field_name == 'Pressure':
        back_path = '../data/converted_data/background_pressure.npz'
        obs_path = '../data/converted_data/obs_pressure.npz'
    elif field_name == 'Velocity':
        back_path = '../data/converted_data/background_velocity.npz'
        obs_path = '../data/converted_data/obs_velocity.npz'


    np.savez_compressed(back_path, u=uTot)
    np.savez_compressed(obs_path, y=y)


def get_positions(filepath):
    ug = vtktools.vtu(filepath + 'LSBU_0.vtu')
    pos = np.array(ug.GetLocations())
    np.savez_compressed('../data/converted_data/positions.npz', pos=pos)

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
    convert(args.filepath, 'Velocity')
    convert(args.filepath, 'Pressure')
    get_positions(args.filepath)
