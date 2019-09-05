#!/usr/bin/env python

import numpy as np
import os
import vtktools
import argparse
from scipy import linalg as LA


def convert(filepath, field_name='Tracer', ntime=988):
    # ntime = 988  # timesteps (there are actually 989 but divide by 2 to split data evenly)
    features = 100040
    uTot = np.zeros([ntime // 2, features])
    y = np.zeros([ntime // 2, features])
    for i in range(ntime):
        print("Processing file " + filename)
        try:
            filename = filepath + 'LSBU_' + str(i) + '.vtu'
        except:
            print(filename, " does not exist or is not a valid vtu file")
            return
        ug = vtktools.vtu(filename)
        if i < ntime / 2:
            if field_name != 'Velocity':
                ui = ug.GetScalarField(field_name)

                # Create the 2D u matrix, [ntime,n]
                # uTot = np.vstack([uTot, ui]) if uTot.size else ui # model
                uTot[i % (ntime // 2)] = ui
            else:
                ui = np.array(ug.GetVectorField(field_name))
                ui = LA.norm(ui, 2, axis=1)

                uTot[i % (ntime // 2)] = ui
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
    print(uTot)
    print(y)

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
    ug = vtktools.vtu(filepath + 'LSBU_0.vtu')
    ug.AddScalarField('x', np.mean(uTot, axis=0))
    ug.AddScalarField('y', np.mean(y, axis=0))
    results_filepath = filepath + 'LSBU_0_results' + field_name + '.vtu'
    if not os.path.exists(results_filepath):
        ug.Write(filepath + 'LSBU_0_results' + field_name + '.vtu')


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
    if args.vel:
        convert(args.filepath, 'Velocity', ntime=args.ntime)
    # convert(args.filepath, 'Pressure')
    if args.pos:
        get_positions(args.filepath)
