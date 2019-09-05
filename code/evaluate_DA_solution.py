import numpy as np
import vtktools
import os
from numpy import linalg as LA

# import vtktools
import argparse


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
    parser.add_argument('-results',
                        '--results_filepath',
                        help='provide xDA solution filepath')
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parser()
    if os.path.basename(args.results_filepath)[-3:] == "vtu":
        try:
            ug = vtktools.vtu(args.results)
        except:
            print("invalid vtu file")
        xDA = ug.GetScalarField('xDA')
        xB = ug.GetScalarField('xB')
        y = ug.GetScalarField('y')
        evaluate_DA_solution(xDA, xB, y)
    elif os.path.basename(args.results_filepath)[-3:] == "npz":
        result = np.load(args.results_filepath)
        print('L2 norm of the error in DA solution', result['control'])
        print('L2 norm of the error in DA solution', result['result'])
        print('Time taken to run DA', result['time'])
        xDA = result['xDA']
        y = np.load(args.y_filepath)['y']
        xB = np.load(args.xB_filepath)['x']
        evaluate_DA_solution(xDA, xB, y)
    else:
        print("Invalid File")
