from BasicCodes.parser import parsing
from BasicCodes.visualization import *
import math

import numpy as np


def getHPWL(net_list, cell_list, pad_list):
    # Enter your code here to get the HPWL
    HPWL = 0
    # get the coordinates
    cell_coordinates = []
    pad_coordinates = []

    for i in range(len(cell_list)):
        theCell: data_structures.Cell = cell_list[i]
        cell_coordinate = [int(theCell.coordinate[0] * 10), int(theCell.coordinate[1] * 10)]
        cell_coordinates.append(cell_coordinate)
    for i in range(len(pad_list)):
        thePad: data_structures.Pad = pad_list[i]
        pad_coordinate = [int(thePad.coordinate[0] * 10), int(thePad.coordinate[1] * 10)]
        pad_coordinates.append(pad_coordinate)


    #################################
    # net로 서로 연결된 pad 또는 cell의 정보를 받아옴.
    cell_coordinates_for_nets = []
    pad_coordinates_for_nets = []

    for i in range(len(net_list)):
        net_i: data_structures.Net = net_list[i]
        connected_cell_num = len(net_i.connected_cells)
        cell_coordinates_for_net_i = []
        for j in range(connected_cell_num):
            theCell: data_structures.Cell = cell_list[net_i.connected_cells[j] - 1]
            cell_coordinate = [int(theCell.coordinate[0] * 10),
                               int(theCell.coordinate[1] * 10)]  # scaled gate coordinate
            cell_coordinates_for_net_i.append(cell_coordinate)
        cell_coordinates_for_nets.append(cell_coordinates_for_net_i)

        connected_pad_num = len(net_i.connected_pad)
        pad_coordinates_for_net_i = []
        for j in range(connected_pad_num):
            thePad: data_structures.Pad = pad_list[net_i.connected_pad[j] - 1]
            pad_coordinate = [int(thePad.coordinate[0] * 10), int(thePad.coordinate[1] * 10)]  # scaled coordinate
            pad_coordinates_for_net_i.append(pad_coordinate)
        pad_coordinates_for_nets.append(pad_coordinates_for_net_i)

    # delete the completely used variables
    del cell_coordinates_for_net_i, cell_coordinate, connected_cell_num, theCell, i, j

    ##################

    X = []
    Y = []

    for i in range(len(net_list)):
        cell_coordinates_for_net_i = cell_coordinates_for_nets[i]
        pad_coordinates_for_net_i = pad_coordinates_for_nets[i]

        for j in range(len(cell_coordinates_for_net_i)):
            X.append(cell_coordinates_for_net_i[j][0])
            Y.append(cell_coordinates_for_net_i[j][1])
        for j in range(len(pad_coordinates_for_net_i)):
            X.append(pad_coordinates_for_net_i[j][0])
            Y.append(pad_coordinates_for_net_i[j][1])
        HPWL += max(X) - min(X) + max(X) - min(Y)
        X = []
        Y = []

    return HPWL


if __name__ == '__main__':
    filename = "benchmarks/toy1"
    net_list, cell_list, pad_list = parsing(filename)

    draw_window(pad_list, cell_list, net_list, ver="visualization")
    HPWL = getHPWL(net_list, cell_list, pad_list)
    print("Initial HPWL:", HPWL)

