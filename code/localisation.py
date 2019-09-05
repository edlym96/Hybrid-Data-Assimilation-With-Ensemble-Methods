#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import argparse
from scipy.sparse.linalg import svds
import os


def localise_h(x_positions, y_positions, cutoff, rh=5):
    # Checks if the full-rank localisation matrix already exists. If not, build it
    if os.path.exists('../data/converted_data/Ch.npy'):
        C = np.load('../data/converted_data/Ch.npy')
    else:
        positions = np.stack((x_positions, y_positions), axis=1)
        print(positions.shape)
        C = np.zeros((positions.shape[0], positions.shape[0]))
        # t = time.time()
        for i in range(positions.shape[0]):
            print("calculating horizontal localisation for ", i)
            # coord = positions[i]
            tmp = np.tile(positions[i], (positions.shape[0] - i, 1))
            s = LA.norm(positions[i:] - tmp, ord=2, axis=1)
            s[s <= cutoff / 2] = 1
            s[s > cutoff] = 0
            s[(s > cutoff / 2) & (s <= cutoff)] = 0.5 * (
                    1 + np.cos((2 * np.pi * (s[(s > cutoff / 2) & (s <= cutoff)] - cutoff / 2)) / cutoff))
            C[i, i:] = s
            C[i:, i] = s
        print("saving matrix Ch...")
        np.save('../data/converted_data/Ch.npy', C)

    # Perform the decomposition on the localisation matrix and take the first rh EOFs
    V, W, _ = svds(C, k=rh)
    ans = np.matmul(V, np.diag(W))
    explained_variance = np.var(ans, axis=0) / np.var(C, axis=0).sum()
    print(explained_variance)
    np.save('../data/results/explained_variance_Ch.npy', explained_variance)
    print(ans.shape)
    return ans


def localise_v(z_positions, scale, rv=5):
    # Checks if the full-rank localisation matrix already exists. If not, build it
    if os.path.exists('../data/converted_data/Cv.npy'):
        C = np.load('../data/converted_data/Cv.npy')
    else:
        C = np.zeros((z_positions.shape[0], z_positions.shape[0]))
        for i in range(z_positions.shape[0]):
            # coord = z_positions[i]
            print("calculating vertical localisation for ", i)
            tmp = np.tile(z_positions[i], (z_positions.shape[0] - i))
            s = np.abs(z_positions[i:] - tmp)
            s = 1 / (1 + (s / scale) ** 2)
            C[i, i:] = s
            C[i:, i] = s
        print("saving matrix Cv...")
        C = np.save('../data/converted_data/Cv.npy', C)

        # Perform the decomposition on the localisation matrix and take the first rv EOFs
    V, W, _ = svds(C, k=rv)
    ans = np.matmul(V, np.diag(W))
    explained_variance = np.var(ans, axis=0) / np.var(C, axis=0).sum()
    print(explained_variance)
    np.save('../data/results/explained_variance_Cv.npy', explained_variance)
    return ans


def arg_parser():
    parser = argparse.ArgumentParser(description='Calculate DA Solution')
    parser.add_argument('-posp',
                        '--pos_filepath',
                        default="../data/converted_data/positions.npz",
                        help='provide file path for position matrix'
                        )
    parser.add_argument('-rh',
                        default=3,
                        help="Horizontal EOF modes"
                        )
    parser.add_argument('-rv',
                        default=3,
                        help='Vertical EOF modes'
                        )
    parser.add_argument('-h_cutoff',
                        default=200,
                        help='Horizontal cutoff'
                        )

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()
    try:
        x_pos = np.load(args.pos_filepath)['pos'][:, 0]
        y_pos = np.load(args.pos_filepath)['pos'][:, 1]
        z_pos = np.load(args.pos_filepath)['pos'][:, 2]
    except:
        print("Position filepath not found. Please run convert_vtu.py with pos flag on")

    C = localise_h(x_pos, y_pos, int(args.h_cutoff), int(args.rh))  # Ch is (100000 , 50)
    np.savez_compressed('../data/converted_data/reduced_localisation_h'+str(args.rh)+'_'+str(args.h_cutoff)+'.npz',C=C)
    C = None
    C = localise_v(z_pos, 10, int(args.rv))
    np.savez_compressed('../data/converted_data/reduced_localisation_v'+str(args.rv)+'.npz', C=C)
