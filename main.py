from BasicCodes.parser import parsing
from BasicCodes.visualization import *
import math

import numpy as np


def getHPWL(net_list, cell_list, pad_list):
    # Enter your code here to get the HPWL
    X = []
    Y = []
    for i in range(len(cell_list)):
        X.append(cell_list[i].coordinate[0])
        Y.append(cell_list[i].coordinate[1])
    HPWL = max(X) - min(X) + max(Y) - min(Y)

    return HPWL


if __name__ == '__main__':
    filename = "benchmarks/toy1"
    net_list, cell_list, pad_list = parsing(filename)

    draw_window(pad_list, cell_list, net_list, ver="visualization")
    HPWL = getHPWL(net_list, cell_list, pad_list)
    print("Initial HPWL:", HPWL)

