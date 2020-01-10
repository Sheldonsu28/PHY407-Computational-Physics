import struct
from copy import deepcopy
import numpy as np

def Question_1b(filename):
    """
    This is the implementation for question 1b
    """
    map_array = np.zeros((1201, 1201))
    file = open(filename, 'rb')
    for y in range(0, 1201):                                    # Read the file into the data array
        for x in range(0, 1201):
            data = file.read(2)
            a = struct.unpack('>h', data)[0]
            if a >= 0:
                map_array[y, x] = a
            else:
                map_array[y, x] = 0
    file.close()
    gradient = calculate_gradient(deepcopy(map_array))               # Calculate the I value suing the gradient

    return gradient, map_array






def calculate_gradient(map_data):
    """
    This function calculates the light intensity for each of the w(x, y)
    """
    dx = []
    dy = []
    for x in range(0, 1201, 10):
        for y in range(0, 1201, 10):

            if x == 0 or x > 1199:               # If x is on the edge, use forward/backwrd difference
                if x == 0:
                    partial_x = derivitive_forward(map_data[x][y], map_data[x + 1][y], 420)

                else:
                    partial_x = derivitive_backward(map_data[x][y], map_data[x - 1][y], 420)

            else:
                partial_x = derivitive_central(map_data[x - 1][y], map_data[x + 1][y], 420)

            if y == 0 or y > 1199:                # If y is on the edge, use forward/backward difference
                if y == 0:
                    partial_y = derivitive_forward(map_data[x][y], map_data[x][y + 1], 420)
                else:
                    partial_y = derivitive_backward(map_data[x][y], map_data[x][y - 1], 420)
            else:
                partial_y = derivitive_central(map_data[x][y - 1], map_data[x][y + 1], 420)

            dx.append(partial_x)
            dy.append(partial_y)

    return (dx, dy)



def derivitive_forward(i,i_1,delta):
    """
    This function implements the forward difference method
    """
    return (i_1 - i)/delta

def derivitive_backward(i,i_1,delta):
    """
    This function implements the backward difference method
    """
    return (i - i_1)/delta

def derivitive_central(i_minus,i_plus,delta):
    """
    This function implements the central difference method
    """
    return (i_plus - i_minus) / (2 * delta)
