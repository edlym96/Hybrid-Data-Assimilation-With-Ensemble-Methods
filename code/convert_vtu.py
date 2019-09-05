#!/usr/bin/env python

import numpy as np
import os
import vtktools
import argparse
from scipy import linalg as LA


def convert(filepath, field_name='Tracer', ntime=988):
    # ntime = 988  # timesteps (there are actually 989 but divide by 2 to split data evenly)
    features = 100040
    xTot = np.zeros([ntime // 2, features])
    y = np.zeros([ntime // 2, features])
    for i in range(ntime):
        try:
            filename = filepath + 'LSBU_' + str(i) + '.vtu'
        except:
            print(filename, " does not exist or is not a valid vtu file")
            return
        print("Processing file " + filename)
        ug = vtktools.vtu(filename)
        if i < ntime / 2:
            if field_name != 'Velocity':
                xi = ug.GetScalarField(field_name)

                # Create the 2D u matrix, [ntime,n]
                # xTot = np.vstack([xTot, xi]) if xTot.size else xi # model
                xTot[i % (ntime // 2)] = xi
            else:
                xi = np.array(ug.GetVectorField(field_name))
                xi = LA.norm(xi, 2, axis=1)

                xTot[i % (ntime // 2)] = xi
        else:
            if field_name != 'Velocity':
                yi = ug.GetScalarField(field_name)
                # y = np.vstack([y,yi]) if y.size else yi # observations
                y[i % (ntime // 2)] = yi
            else:
                yi = ug.GetVectorField(field_name)
                yi = LA.norm(yi, 2, axis=1)

                y[i % (ntime // 2)] = yi

    print(ug.GetFieldNames())
    print(xTot.shape)
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

    np.savez_compressed(back_path, x=xTot)
    np.savez_compressed(obs_path, y=y)


def get_positions(filepath):
    try:
        ug = vtktools.vtu(filepath + 'LSBU_0.vtu')
    except:
        print("LSBU_0.vtu does not exist or is not a valid vtu file")
        return
    pos = np.array(ug.GetLocations())
    np.savez_compressed('../data/converted_data/positions.npz', pos=pos)


def arg_parser():
    parser = argparse.ArgumentParser(description='Convert vtu to npz')
    parser.add_argument('-fp',
                        '--filepath',
                        default="../data/small3DLSBU/",
                        help='provide file path data')
    parser.add_argument('-ntime',
                        default=988)
    parser.add_argument('-pos',
                        '--positions',
                        action='store_true')
    parser.add_argument('-vel',
                        '--velocity',
                        action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()

    convert(args.filepath, ntime=args.ntime)
    if args.velocity:
        convert(args.filepath, 'Velocity', ntime=args.ntime)
    # convert(args.filepath, 'Pressure')
    if args.positions:
        get_positions(args.filepath)
