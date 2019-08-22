#!/usr/bin/env python

import numpy as np
import scipy.linalg as LA
import argparse

def localise_h(x_positions, y_positions, cutoff, rh=5):
    # Lh = max(positions) - min(positions)
    # C = np.zeros([len(positions), len(positions)])
    # TODO: MOVE CONSTRUCTION OF LOCALISATION MATRIX TO OPTIMAL COVARIANCE
    # downsampled_x = x_positions[::scale]
    # downsampled_y = y_positions[::scale]

    # positions = np.stack((downsampled_x,downsampled_y), axis=1)
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
        """for j in range(i,positions.shape[0]):
			print(i,j)
			s = LA.norm(positions[j]-coord, ord=2)
			if s <= cutoff / 2:
				C[i,j] = C[j,i]= 1
			elif s > cutoff:
				C[i,j] = C[j,i] = 0
			else:
				C[i,j] = C[j,i] = 0.5*(1+math.cos((2*math.pi*(s-cutoff/2))/cutoff))
		"""
    W, V = LA.eigh(C, eigvals=(C.shape[0] - rh, C.shape[0] - 1))
    print(C)
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
    ans = np.matmul(V, np.diag(W))
    return ans


def localise_v(z_positions, scale, rv=5):
    C = np.zeros((z_positions.shape[0], z_positions.shape[0]))
    for i in range(z_positions.shape[0]):
        # coord = z_positions[i]
        print("calculating vertical localisation for ", i)
        tmp = np.tile(z_positions[i], (z_positions.shape[0] - i, 1))
        s = np.abs(z_positions[i:] - tmp)
        s = 1 / (1 + (s / scale) ** 2)
        C[i, i:] = s
        C[i:, i] = s
        """for j in range(i,z_positions.shape[0]):
            dz = abs(z_positions[j]-coord)
            C[i, j] = C[j, i] = 1/(1+(dz/scale)**2)
        """
    W, V = LA.eigh(C, eigvals=(C.shape[0] - rv, C.shape[0] - 1))
    idx = W.argsort()[::-1]
    # W = np.sqrt(W[idx])
    W = W[idx]
    V = V[:, idx]
    return np.matmul(V, np.diag(W))


def arg_parser():
    parser = argparse.ArgumentParser(description='Calculate DA Solution')
    parser.add_argument('-posp',
                        '--pos_filepath',
                        default="../data/converted_data/positions.npz",
                        help='provide file path for position matrix'
                        )
    parser.add_argument('-rj',
                        default=5,
                        help="Horizontal EOF modes"
                        )
    parser.add_argument('-rv',
                        default=5,
                        help='Vertical EOF modes'
                        )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()
    x_pos = np.load(args.pos_filepath)['pos'][:, 0]
    y_pos = np.load(args.pos_filepath)['pos'][:, 1]
    z_pos = np.load(args.pos_filepath)['pos'][:, 2]
    C = localise_h(x_pos, y_pos, 200)  # Ch is (100000 , 50)
    np.savez_compressed('../data/converted_data/localisation_h.npz',C=C)
    C = None
    C = localise_v(z_pos, 10)
    np.savez_compressed('../data/converted_data/localisation_v.npz', C=C)
