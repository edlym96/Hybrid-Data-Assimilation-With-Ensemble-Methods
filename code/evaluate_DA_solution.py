import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
import os

from numpy import linalg as LA

import math

import sys
#import vtktools
import argparse
#sys.path.append('fluidity-master')


def evaluate_DA_solution(xDA, xB, y):
    errxB = y - xB
    MSExb = LA.norm(errxB, 2) / LA.norm(y, 2)
    print('L2 norm of the background error', MSExb, '\n')

    errxDA = y - xDA
    MSExDA = LA.norm(errxDA, 2) / LA.norm(y, 2)
    print('L2 norm of the error in DA solution', MSExDA, '\n')

    return MSExDA

def arg_parser():
    parser = argparse.ArgumentParser(description='Evaluate DA solution')
    parser.add_argument('-xDA',
                        '--xDA_filepath',
                        help='provide xDA solution filepath')
    parser.add_argument('-xB',
                        '--xB_filepath',
                        default="../data/converted_data/background_state.npz",
                        help='provide file path for background data')
    parser.add_argument('-y',
                        '--y_filepath',
                        default="../data/converted_data/observations.npz",
                        help='provide file path for observation data')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()

    evaluate_DA_solution(args.xDa, args.xB, args.y)