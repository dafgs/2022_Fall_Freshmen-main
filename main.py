from BasicCodes.parser import parsing
from BasicCodes.visualization import *
import math

import numpy as np


def getHPWL(net_list, cell_list, pad_list):
    # Enter your code here to get the HPWL
    net_coordinates = np.zeros((len(net_list), 2), dtype=np.float64)
    X = []
    Y = []
    for i in range(len(net_list)):
        net_coordinates[i, 0] = net_list[i].coordinate[0]
        X.append(net_list[i].coordinate[0])
        net_coordinates[i, 1] = net_list[i].coordinate[1]
        Y.append(net_list[i].coordinate[1])

    HPWL = max(X) - min(X) + max(Y) - min(Y)

    '''
    pad_coordinates = np.zeros((len(net_list), 2), dtype=np.float64)
    X = []
    Y = []
    for i in range(len(pad_list)):
        pad_coordinates[i, 0] = pad_list[i].coordinate[0]
        X.append(pad_list[i].coordinate[0])
        pad_coordinates[i, 1] = pad_list[i].coordinate[1]
        Y.append(pad_list[i].coordinate[1])

    HPWL = max(X) - min(X) + max(Y) - min(Y)
    '''

    return HPWL


if __name__ == '__main__':
    filename = "benchmarks/toy1"
    net_list, cell_list, pad_list = parsing(filename)

    draw_window(pad_list, cell_list, net_list, ver="visualization")
    HPWL = getHPWL(net_list, cell_list, pad_list)
    print("Initial HPWL:", HPWL)

