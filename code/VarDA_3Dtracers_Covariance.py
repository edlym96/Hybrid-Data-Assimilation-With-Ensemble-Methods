#!/usr/bin/env python3
import numpy as np
from scipy.optimize import minimize
import time
import os

from scipy import linalg as LA

import vtktools
import argparse


# sys.path.append('fluidity-master')


def build_DA_solution(xB_filepath, y_filepath, V_filepath, ntime=989 // 2, h_localisation=None,
                      v_localisation=None):
    try:
        xB = np.transpose(
            np.load(xB_filepath)['x'])  # Transpose from 989x100040 to 100040x989
        print("xB", xB.shape)
    except:
        print("Background error covariance matrix not found. Please run convert_vtu.py")
        return

    n = xB.shape[0]
    try:
        y = np.transpose(
            np.load(y_filepath)['y'])  # Transpose from 989x100040 to 100040x989
        print("y", y.shape)
    except:
        print("Observation matrix not found. Please run convert_vtu.py")
        return

    try:
        V = np.load(V_filepath)['V']
    except:
        print("V matrix not found. Please run optimal_covariance.py")
        return

    lam = 0.1e-60
    R = lam * 0.9
    x0 = np.zeros(n)

    # Background state in obervation space. H is taken as identity matrix so a direct copy can be used.
    HxB = xB.copy()

    # Misfit calculated by subtracting
    d = np.subtract(y, HxB)

    if h_localisation and v_localisation:
        try:
            Ch = np.load(h_localisation)['C']
            Cv = np.load(v_localisation)['C']
        except:
            print("Localisation matrix not found. Please run localisation.py")
            return
        V_new = np.zeros([V.shape[0], V.shape[1] * Ch.shape[1] * Cv.shape[1]])
        for i in range(V.shape[1]):
            tmp = np.tile(V[:, i], (Ch.shape[1], 1)).transpose()  # Shape (100000,50)
            tmp = np.multiply(Ch, tmp)
            for j in range(tmp.shape[1]):
                tmp2 = np.tile(tmp[:, j], (Cv.shape[1], 1)).transpose()  # Shape (100000,40)
                tmp2 = np.multiply(Cv, tmp2)
                V_new[:, (i * Ch.shape[1] + j) * Cv.shape[1]:(i * Ch.shape[1] + j + 1) * Cv.shape[1]] = tmp2
        V = V_new
    elif h_localisation and not v_localisation:
        try:
            Ch = np.load(h_localisation)['C']
        except:
            print("Horizontal Localisation matrix not found. Please run localisation.py")
            return
        V_new = np.zeros([V.shape[0], V.shape[1] * Ch.shape[1]])
        for i in range(V.shape[1]):
            tmp = np.tile(V[:, i], (Ch.shape[1], 1)).transpose()  # Shape (100000,50)
            tmp = np.multiply(Ch, tmp)
            V_new[:, i * Ch.shape[1]:(i + 1) * Ch.shape[1]] = tmp
        V = V_new
    elif v_localisation and not h_localisation:
        try:
            Cv = np.load(v_localisation)['C']
        except:
            print("Vertical Localisation matrix not found. Please run localisation.py")
            return
        V_new = np.zeros([V.shape[0], V.shape[1] * Cv.shape[1]])
        for i in range(V.shape[1]):
            tmp = np.tile(V[:, i], (Cv.shape[1], 1)).transpose()  # Shape (100000,50)
            tmp = np.multiply(Cv, tmp)
            V_new[:, i * Cv.shape[1]:(i + 1) * Cv.shape[1]] = tmp
        V = V_new

    Vin = np.linalg.pinv(V)

    print("Vin", Vin.shape)

    v0 = np.dot(Vin,
                x0)  # Take x0 from physical to reduced space by dotting with inverse of reduced V
    print("v0 ", v0.shape)
    VT = np.transpose(V)

    deltaxDA = np.zeros((n, ntime))

    t = time.time()
    for i in range(ntime):
        print("Processing timestep no: ", i)

        di = d[:, i].reshape([n, 1])

        # Cost function J
        def J(v):
            # v = v.reshape((V.shape[1],ntime)) # need to reshape because optimize flattens input
            v = v.reshape((V.shape[1], 1))
            vT = np.transpose(v)
            vTv = np.dot(vT, v)
            Vv = np.dot(V, v).reshape([n, 1])
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
            Vv = np.dot(V, v).reshape([n, 1])
            Jmis = np.subtract(Vv, di)
            invR = 1 / R
            gg2 = np.multiply(invR, np.dot(VT, Jmis))  # VT[501x100040] Jmis[989x100040]
            ggJ = v + gg2
            # return ggJ.flatten()
            return ggJ

        res = minimize(J, v0, method='L-BFGS-B', jac=gradJ, options={'disp': False})
        vDA = np.reshape(res.x, (V.shape[1], 1))
        deltaxDA[:, i] = np.dot(V, vDA).flatten()  # take vDA from the reduced space back to x-space
        v0 = vDA  # FGAT step

    elapsed = time.time() - t
    print('elapsed', elapsed, 'seconds\n')
    xDA = xB + deltaxDA

    MSExb, MSExDA = evaluate_DA_solution(xDA, xB, y)

    results_filename = os.path.basename(V_filepath)

    save_DA_solution(xDA, deltaxDA, y, MSExDA, MSExb, results_filename, h_localisation, v_localisation, elapsed)


def save_DA_solution(xDA, deltaxDA, y, MSE, MSExb, filename, h_localisation, v_localisation, elapsed):
    if not os.path.exists('../data/results'):
        os.makedirs('../data/results')
    print("Saving results to " + filename + "...")
    try:
        ug = vtktools.vtu('../data/small3DLSBU/LSBU_0.vtu')
    except:
        print("Results matrix to record results not found. Please run convert_vtu.py again")
        return
    ug.AddScalarField('abs y-xDA', np.mean(np.abs(y - xDA), axis=1))
    ug.AddScalarField('abs y-xB', np.mean(np.abs(y - xDA + deltaxDA), axis=1))
    ug.AddScalarField('deltaxDA', np.mean(deltaxDA, axis=1))
    ug.AddScalarField('y-xDA', np.mean(y - xDA, axis=1))
    ug.AddScalarField('y-xB', np.mean(y - xDA + deltaxDA, axis=1))
    ug.AddScalarField('xDA', np.mean(xDA, axis=1))
    ug.AddScalarField('y', np.mean(y, axis=1))
    ug.AddScalarField('xB', np.mean(xDA - deltaxDA, axis=1))
    path = "../data/results/Results"
    if h_localisation or v_localisation:
        # path += "Localisation"
        path += "LocalisationFGAT"
        if h_localisation:
            path += h_localisation.replace("../data/converted_data/reduced_localisation_h", 'rh').replace('.npz', '')
        else:
            path += 'rh0'
        if v_localisation:
            # path += 'rv'+str(rv)
            path += v_localisation.replace("../data/converted_data/reduced_localisation_v", 'rv').replace('.npz', '')
        else:
            path += 'rv0'
    path += filename
    print("Saving results to " + path + "...")
    np.savez_compressed(path, result=MSE, time=elapsed, control=MSExb, xDA=xDA)  # SAve MSE, control and time taken
    ug.Write(path.replace('.npz', '.vtu'))  # Save the results of DA to vtu


def evaluate_DA_solution(xDA, xB, y):
    errxB = y - xB
    MSExb = LA.norm(errxB, 2) / LA.norm(y, 2)
    print('L2 norm of the background error', MSExb, '\n')

    errxDA = y - xDA
    MSExDA = LA.norm(errxDA, 2) / LA.norm(y, 2)
    print('L2 norm of the error in DA solution', MSExDA, '\n')

    return MSExb, MSExDA


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
    parser.add_argument('--h_localisation',
                        '-hlocal',
                        default='',
                        help='Horizontal Localisation Matrix Filepath'
                        )
    parser.add_argument('--v_localisation',
                        '-vlocal',
                        default='',
                        help='Vertical Localisation Matrix Filepath'
                        )

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()

    build_DA_solution(args.xB_filepath, args.y_filepath, args.V_filepath, args.ntime,
                      args.h_localisation, args.v_localisation)
