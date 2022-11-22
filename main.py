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
###############################################################
# ASSIGNMENT 2

def solve_AX_b(A,b):
    inverse_A = np.linalg.inv(A)
    #print(inverse_A)
    result = np.dot(b,inverse_A)
    #print(b)
    return result

def QuadraticPlacement(net_list, cell_list, pad_list):

    cell_number = len(cell_list)
    pad_number = len(pad_list)
    # cell 또는 pad가 연결되어 있으면 1 아니면 0인 표. 즉, 가중치 표
    weight = [[0 for _ in range(cell_number+pad_number)] for _ in range(cell_number+pad_number)]

    #################################
    # STEP1: net로 서로 연결된 pad와 cell을 조사하여 가중치 표를 만듦.

    for i in range(len(net_list)):
        net_i: data_structures.Net = net_list[i]

        connected_cell_num = len(net_i.connected_cells)
        connected_pad_num = len(net_i.connected_pad)
        cell_id_and_pad_id_for_net_i = []
        for j in range(connected_cell_num):
            theCell: data_structures.Cell = cell_list[net_i.connected_cells[j] - 1]
            cell_id = int(theCell.id)
            cell_id_and_pad_id_for_net_i.append(cell_id)

        for j in range(connected_pad_num):
            thePad: data_structures.Pad = pad_list[net_i.connected_pad[j] - 1]
            pad_id = int(thePad.pad_id)
            cell_id_and_pad_id_for_net_i.append(pad_id)
        for j in range(len(cell_id_and_pad_id_for_net_i)):
            for k in range(1, len(cell_id_and_pad_id_for_net_i) - j):
                if j<connected_cell_num and j+k<connected_cell_num:
                    weight[cell_id_and_pad_id_for_net_i[j]][cell_id_and_pad_id_for_net_i[j+k]] += 1  # 한 net로 연결되어 있는 Cell들을 콕콕 집어가며 표에 1을 더해줌.
                    weight[cell_id_and_pad_id_for_net_i[j+k]][cell_id_and_pad_id_for_net_i[j]] += 1
                elif j<connected_cell_num and j+k>=connected_cell_num:
                    weight[cell_id_and_pad_id_for_net_i[j]][cell_number+cell_id_and_pad_id_for_net_i[j+k]] += 1  # 한 net로 연결되어 있는 Cell들을 콕콕 집어가며 표에 1을 더해줌.
                    weight[cell_number+cell_id_and_pad_id_for_net_i[j+k]][cell_id_and_pad_id_for_net_i[j]] += 1
                elif j >= connected_cell_num and j+k>connected_cell_num:
                    weight[cell_number+cell_id_and_pad_id_for_net_i[j]][cell_number + cell_id_and_pad_id_for_net_i[j + k]] += 1  # 한 net로 연결되어 있는 Cell들을 콕콕 집어가며 표에 1을 더해줌.
                    weight[cell_number + cell_id_and_pad_id_for_net_i[j + k]][cell_number+cell_id_and_pad_id_for_net_i[j]] += 1

    # delete the completely used variables
    del connected_cell_num, theCell, i, j, k

    # STEP2: 미분을 했을 때 0되게 하는 X를 구하기 위해 행렬 식 만듦
    diagonal = [[0 for _ in range(cell_number)] for _ in range(cell_number)]
    for i in range(cell_number):
        sum = 0
        for j in range(cell_number+pad_number):
            sum += weight[i][j]
        diagonal[i][i] += sum

    Connectivitiy_matrix = [[0 for _ in range(cell_number)] for _ in range(cell_number)]
    for i in range(cell_number):
        for j in range(cell_number):
            Connectivitiy_matrix[i][j] += weight[i][j]

    A = np.matrix(diagonal) - np.matrix(Connectivitiy_matrix)

    # for x,
    b_x = [[0 for _ in range(1)] for _ in range(cell_number)]

    k = 0
    for i in range(pad_number):
        for j in range(cell_number):
            if weight[i+cell_number][j] != 0:
                b_x[k][0] += weight[i+cell_number][j]*pad_list[i].coordinate[0]
                k += 1
                break

    # for y,
    b_y = [[0 for _ in range(1)] for _ in range(cell_number)]

    k = 0
    for i in range(pad_number):
        for j in range(cell_number):
            if weight[i + cell_number][j] != 0:
                b_y[k][0] += weight[i + cell_number][j] * pad_list[i].coordinate[1]
                k += 1
                break

    b_x = np.matrix(b_x)
    b_y = np.matrix(b_y)

    b_x = np.transpose(b_x)
    b_y = np.transpose(b_y)

    # STEP3: 최적의 (x,y)를 구함.
    cells_x = solve_AX_b(A,b_x) # Quadratic Placement의 cell의 x좌표들(세로로 배열되어 있음)
    cells_y = solve_AX_b(A,b_y) # Quadratic Placement의 cell의 y좌표들(세로로 배열되어 있음)

    # STEP4: 위에서 구한 x,y를 바탕으로 cell의 좌표를 바꿔줄것임.
    for i in range(len(cell_list)):
        temp = int(cell_list[i].coordinate[0])
        cell_list[i].coordinate[0] -= temp
        cell_list[i].coordinate[0] += cells_x[0,i]

        temp = int(cell_list[i].coordinate[1])
        cell_list[i].coordinate[1] -= temp
        cell_list[i].coordinate[1] += cells_y[0,i]


if __name__ == '__main__':
    #for filename in ["benchmarks/biomed","benchmarks/industry1","benchmarks/industry2", "benchmarks/primary","benchmarks/struct","benchmarks/toy1", "benchmarks/toy2"]:
        filename = "benchmarks/toy1"
        net_list, cell_list, pad_list = parsing(filename)

        # random placement
        draw_window(pad_list, cell_list, net_list, ver="visualization")
        HPWL = getHPWL(net_list, cell_list, pad_list)
        print("Random Placement HPWL=", HPWL)

        # quadraticPlacement
        QuadraticPlacement(net_list, cell_list, pad_list)
        draw_window(pad_list, cell_list, net_list, ver="visualization")
        HPWL = getHPWL(net_list, cell_list, pad_list)
        print("Quadratic Placement HPWL=", HPWL)

