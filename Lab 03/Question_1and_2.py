import struct
import pylab as plt
from copy import deepcopy
from guassxw import gaussxw
import math as m
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
            map_array[y, x] = struct.unpack('>h', data)[0]
    file.close()
    I = calculate_Intensity(deepcopy(map_array))               # Calculate the I value suing the gradient
    plt.figure(0)
    plt.imshow(map_array, vmin=0, vmax=4400)
    plt.colorbar()
    plt.title("Map data readed from the file")
    plt.legend()
    plt.savefig("1b_picture_Intesity.png")
    plt.show()





def calculate_Intensity(map_data):
    """
    This function calculates the light intensity for each of the w(x, y)
    """
    I = np.empty((1021, 1021))
    for x in range(0, 1021):
        for y in range(0, 1021):

            if x == 0 or x > 1019:               # If x is on the edge, use forward/backwrd difference
                if x == 0:
                    partial_x = derivitive_forward(map_data[x][y], map_data[x + 1][y], 420)
                else:
                    partial_x = derivitive_backward(map_data[x][y], map_data[x - 1][y], 420)
            else:
                partial_x = derivitive_central(map_data[x - 1][y], map_data[x + 1][y], 420)

            if y == 0 or y > 1019:                # If y is on the edge, use forward/backward difference
                if y == 0:
                    partial_y = derivitive_forward(map_data[x][y], map_data[x][y + 1], 420)
                else:
                    partial_y = derivitive_backward(map_data[x][y], map_data[x][y - 1], 420)
            else:
                partial_y = derivitive_central(map_data[x][y - 1], map_data[x][y + 1], 420)

            I[x][y] =-(m.cos(m.pi) * partial_y + m.sin(m.pi) * partial_x)\
                      / m.sqrt((partial_x**2) + (partial_y**2) + 1)

    return I



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


def Question_2a():
    """
    This is the implimentation to question 2a
    """
    wg_8 = []                                                       # Initialize all the parameters
    wg_16 = []
    x_8, w_8 = gaussxw(8)                                           # Data points for N = 8
    x_16, w_16 = gaussxw(16)                                        # Data points for N = 16

    x_partition_8 = 0.5 * 0.01 * x_8 + 0.5 * 0.01
    x_partition_16 = 0.5 * 0.01 * x_16 + 0.5 * 0.01                 # Map the datapoints to domain of the function

    w_partition_8 = 0.5 * 0.01 * w_8
    w_partition_16 = 0.5 * 0.01 * w_16

    for n in range(0, 8):
        wg_8.append(4 * w_partition_8[n] * g(x_partition_8[n], 0.01))
    for r in range(0, 16):
        wg_16.append(4 * w_partition_16[r] * g(x_partition_16[r], 0.01))
    exact_period = 2 * m.pi * m.sqrt(1/12)
    print((sum(wg_16) - exact_period) /exact_period )
    print((sum(wg_8) - exact_period)/exact_period)
    plt.plot(x_partition_8, wg_8, "ro", label="Approximation with N = 8")
    plt.plot(x_partition_16, wg_16, "bo", label="Approximation with N = 16")
    plt.xlabel("Value of x")
    plt.ylabel("Value of y")
    plt.title("Reimann Sum for ntegration")
    plt.savefig("Question2a.png")
    plt.show()


def Question_2c_1():
    exact_period = 2 * m.pi * m.sqrt(1 / 12)                                # Initialize all the parameters
    wg_200 = []                                                             # Partition for N = 200
    wg_100 = []                                                             # Partition for N = 100
    x_200, w_200 = gaussxw(200)                                             # Calculate the datapoints
    x_100, w_100 = gaussxw(100)

    x_partition_200 = 0.5 * 0.01 * x_200 + 0.5 * 0.01
    w_partition_200 = 0.5 * 0.01 * w_200

    x_partition_100 = 0.5 * 0.01 * x_100 + 0.5 * 0.01
    w_partition_100 = 0.5 * 0.01 * w_100

    for n in range(0, 200):                                                # put the sum element into the corresponding list
        wg_200.append(4 * w_partition_200[n] * g(x_partition_200[n], 0.01))

    for b in range(0, 100):
        wg_100.append(4 * w_partition_100[b] * g(x_partition_100[b], 0.01))

    print(abs(sum(wg_200) - exact_period)/exact_period)                    # Calculate the sum
    print((sum(wg_200) - sum(wg_100))/exact_period)
    plt.plot(x_partition_200, wg_200, "ro")
    plt.show()

def Question_2c_2():
    """
    code for question 2c that calculate and plot the graph for question 2c
    """
    x_value = []                                                # Initialize all the parameters
    relativistic_limit = []
    classical_limit = []
    relativistic_period = []
    for r in range(1, 866025404, 10000):                        # Caclulate peroid T for 1 < x < 866025404, step 10000

        x, w = gaussxw(200)
        x = np.array(x)
        w = np.array(w)

        x_partition = 0.5 * r * x + 0.5 * r                     # Similar to previous codes

        w_partition = 0.5 * r * w
        w_partition = np.array(w_partition)

        vectfun = np.vectorize(g, otypes=[np.float], cache=False)  # Optimize using numpy.vectorization

        total = np.sum(4 * w_partition * vectfun(x_partition, r))
        relativistic_period.append(total)
        x_value.append(r)
        relativistic_limit.append(4 * r / 3.0e8)            # Calculate the relativistic limit
        classical_limit.append(2 * m.pi * m.sqrt(1/12))            # Calculate the classical limit

    plt.plot(x_value, relativistic_period, "bo", label="Period Time")       # Plot the data
    plt.plot(x_value, classical_limit, "go", label="Classical limit")
    plt.plot(x_value, relativistic_limit, "ro", label="Relativistic limit")
    plt.xlabel("Initial posistion of x, unit: m")
    plt.ylabel("Period time, unit: s")
    plt.title("Period versus Relativistic and classical limit ")
    plt.legend()
    plt.savefig("Question2c_2.png")
    plt.show()

def g(x,x_0):
    """
    This is the function g's implementation in the lab handout
    """
    c = 3.0e8
    mass = 1
    k = 12
    return (c * ((0.5 * k * ((x_0**2) - (x**2)) * (2 * mass * (c ** 2) + 0.5 * k * ((x_0 ** 2) - (x ** 2)))) /
                 ((mass * c ** 2) + 0.5 * k * ((x_0 ** 2) - (x ** 2))) ** 2) ** 0.5) ** -1





if __name__ == "__main__":
    Question_1b("N19W156.hgt")
    #Question_2a()
    #Question_2c_1()
    #Question_2c_2()
