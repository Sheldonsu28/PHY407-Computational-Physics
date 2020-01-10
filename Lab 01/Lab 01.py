
import pylab as plt
import math
import time
import numpy as np
import scipy.optimize
import random
from typing import Tuple, List


def solution_to_1c(x_initial, y_initial, x_v_inital, y_v_initial, time_step, end_time):
    """
    This programs plots the solution to 1c.
    """
    x_coordinate = [0]                                                                                                  # Add the position of the sun
    y_coordinate = [0]                                                                                                  # size of x_coordinate and y_coordinate should be = to size of time
    x_new_initial = x_initial
    x_new_v_initial = x_v_inital
    y_new_initial = y_initial
    y_new_v_initial = y_v_initial

    iteration = int(end_time / time_step)

    for x in range(0, iteration):

        r = math.sqrt((x_new_initial ** 2) + (y_new_initial ** 2))                                                      # Calculate r

        x = x_new_initial + x_new_v_initial * time_step - ((39.5 * 1 * x_new_initial) / (r ** 3)) * (time_step ** 2)    # Calculate the new step value for x
        y = y_new_initial + y_new_v_initial * time_step - ((39.5 * 1 * y_new_initial) / (r ** 3)) * (time_step ** 2)    # Calculate the new step value for y

        x_v = x_new_v_initial - ((39.5 * 1 * x_new_initial) / (r ** 3)) * time_step                                     # Calculate the new value for x velocity
        y_v = y_new_v_initial - ((39.5 * 1 * y_new_initial) / (r ** 3)) * time_step                                     # Calculate the new value for y velocity

        y_coordinate.append(y)                                                                                          # Add new coordinates to list
        x_coordinate.append(x)

        x_new_initial = x                                                                                               # Update starting position values
        y_new_initial = y

        x_new_v_initial = x_v                                                                                           # Update starting velocity values
        y_new_v_initial = y_v

    plt.plot(x_coordinate, y_coordinate, "bo")                                                                          # Plot them to the graph
    plt.xlabel("X position, Unit: AU")
    plt.ylabel("Y position, Unit: AU")
    plt.savefig("1c_ouptput.png")
    plt.show()


def solution_to_1d(x_initial, y_initial, x_v_inital, y_v_initial, time_step, end_time):
    """
    This programs plots the solution to 1d, where the gravitational forces factor predict by general relativity is added
    to the equation.
    """
    x_coordinate = [0]                                                                                                  # Add the position of the sun
    y_coordinate = [0]                                                                                                  # size of x_coordinate and y_coordinate should be = to size of time
    x_new_initial = x_initial
    x_new_v_initial = x_v_inital
    y_new_initial = y_initial
    y_new_v_initial = y_v_initial

    iteration = int(end_time / time_step)

    for x in range(0, iteration):

        r = math.sqrt((x_new_initial ** 2) + (y_new_initial ** 2))                                                      # Calculate r
        rc = (1 + (0.01 / (r**2)))                                                                                      # Calculate general relativity coefficent

        x = x_new_initial + x_new_v_initial * time_step - ((39.5 * 1 * x_new_initial) * rc / (r ** 3)) *(time_step ** 2)# Calculate the new step value for x
        y = y_new_initial + y_new_v_initial * time_step - ((39.5 * 1 * y_new_initial) * rc / (r ** 3)) *(time_step ** 2)# Calculate the new step value for y

        x_v = x_new_v_initial - ((39.5 * 1 * x_new_initial) * rc / (r ** 3)) * time_step                                # Calculate the new value for x velocity
        y_v = y_new_v_initial - ((39.5 * 1 * y_new_initial) * rc / (r ** 3)) * time_step                                # Calculate the new value for y velocity

        y_coordinate.append(y)                                                                                          # Add new coordinates to list
        x_coordinate.append(x)

        x_new_initial = x                                                                                               # Update starting position values
        y_new_initial = y

        x_new_v_initial = x_v                                                                                           # Update starting velocity values
        y_new_v_initial = y_v

    plt.plot(x_coordinate, y_coordinate, "bo")                                                                          # Plot them to the graph
    plt.xlabel("X position, Unit: AU")
    plt.ylabel("Y position, Unit: AU")
    plt.savefig("1d_ouptput.png")
    plt.show()

def solution_to_problem_3_N():
    """
    This function measures the time for matrix muiltiplication on matrix of size N, the function will produce a graph
    with respect to N, apology in advance, because it might take 1 hour or more to finish running because my
    implementation of matrix muiltiplication is a runtime disaster

    """
    # Below is the code is for ploting time as function of N
    N = []
    time_spend_dot = []
    time_spend_multiple = []
    for x in range(2, 216, 1):                  # Matrix size from 2 to 216
        time_average_dot = []
        time_average_multiple = []
        for y in range(0, 5):                  # for each size, we run 5 times to get more consistant data

            A = np.ones([x, x], float) * 5      # generate the matrices
            B = np.ones([x, x], float) * 8

            inital_time = time.time() * 1000            # count initial time for dot product in milisecond

            C = np.dot(A, B)

            final_time = time.time() * 1000             # count final time for for product in milisecond

            inital_time_multiply = time.time()*1000     # count initial time for my own matrix product in milisecond

            D = matrix_muiltiplication(A, B)

            final_time_multiply = time.time()*1000      # count final time for my own matrix product in milisecond

            time_average_dot.append(final_time - inital_time)                       # Find the time spent for all muiltiplications
            time_average_multiple.append(final_time_multiply - inital_time_multiply)

        N.append(x)
        time_spend_dot.append(sum(time_average_dot)/len(time_average_dot))          # Find the averages of the 10 runs
        time_spend_multiple.append((sum(time_average_multiple) / len(time_average_multiple)))

    constants_dot, variant_dot = scipy.optimize.curve_fit(helper_Q3, N, time_spend_dot, (0, 0.0))       # Curve fit for both matrix muiltiplications
    constants_multiple, variant_multiple = scipy.optimize.curve_fit(helper_Q3, N, time_spend_multiple, (0, 0.0))

    time_fit_dot = helper_Q3(N, constants_dot[0], constants_dot[1])                 # calculate the curve using the value calculated
    time_fit_multiple = helper_Q3(N, constants_multiple[0], constants_multiple[1])

    plt.plot(N, time_fit_dot, "bo", label = "My implementation of matrix muiltiplication")
    plt.plot(N, time_fit_multiple, "ro", label = "Numpy.dot")
    plt.xlabel("Matrix size")
    plt.ylabel("Time spent")
    plt.title("Time spent for matrix mutilplication as a function of N")
    plt.legend()
    plt.savefig("problem_3_N.png")
    plt.show()


def solution_to_problem_3_N_cube():
        """
        This function measures the time for matrix muiltiplication on matrix of size N, the function will produce a graph
        with respect to N cube
        """
        # Below is the code is for ploting time as function of N
        n_1 = []
        time_spend_dot_1 = []
        time_spend_multiple_1 = []
        for x in range(2, 6, 1):                                                                        # Matrix size from 2**3 to 6**3
            time_average_dot_1 = []
            time_average_multiple_1 = []
            for y in range(0, 5):                                                                       # for each size, we run 5 times to get more consistant data

                a_1 = np.ones([x**3, x**3], float) * 5                                                  # generate the matrices
                b_1 = np.ones([x**3, x**3], float) * 8

                inital_time_1 = time.time() * 1000                                                      # count initial time for dot product in milisecond

                c_1 = np.dot(a_1, b_1)

                final_time_1 = time.time() * 1000                                                       # count final time for for product in milisecond

                inital_time_multiply_1 = time.time() * 1000                                             # count initial time for my own matrix product in milisecond

                d_1 = matrix_muiltiplication(a_1, b_1)

                final_time_multiply_1 = time.time() * 1000                                              # count final time for my own matrix product in milisecond

                time_average_dot_1.append(final_time_1 - inital_time_1)                                 # Find the time spent for all muiltiplications
                time_average_multiple_1.append(final_time_multiply_1 - inital_time_multiply_1)

            n_1.append(x**3)
            time_spend_dot_1.append(sum(time_average_dot_1) / len(time_average_dot_1))                  # Find the averages of the 10 runs
            time_spend_multiple_1.append((sum(time_average_multiple_1) / len(time_average_multiple_1)))

        constants_dot_1, variant_dot_1 = scipy.optimize.curve_fit(helper_Q3, n_1, time_spend_dot_1, (0, 0.0))  # Curve fit for both matrix muiltiplications
        constants_multiple_1, variant_multiple_1 = scipy.optimize.curve_fit(helper_Q3, n_1, time_spend_multiple_1, (0, 0.0))

        time_fit_dot_1 = helper_Q3(n_1, constants_dot_1[0], constants_dot_1[1])                               # calculate the curve using the value calculated
        time_fit_multiple_1 = helper_Q3(n_1, constants_multiple_1[0], constants_multiple_1[1])

        plt.plot(n_1, time_fit_dot_1, "b", label = "My implementation")
        plt.plot(n_1, time_fit_multiple_1, "r", label = "Numpy.dot")
        plt.xlabel("Matrix size, N**3")
        plt.ylabel("Time spent")
        plt.title("Time spent for matrix mutilplication as a function of N")
        plt.legend()
        plt.savefig("problem_3_n_cube.png")
        plt.show()

def helper_Q3(x, a, b):

    return a * np.power(x, 3) + b

def matrix_muiltiplication(a,b):
    """ My own implementation of matrix muiltiplcation algorithm. Assume the size of the matrices are both n*n
        This is a really bad implementation
    """
    return_matrix = []
    size = len(a)
    for i in range(0, size):        # Create a matrix that the size is n*n (we do all the addition on this matrix)
        return_matrix.append([])

    for columns in return_matrix:   # Fill the matrix with zeros
        for x in range(0, size):
            columns.append(0)
    for x in range(0, size):        # iterate through a's rows
        for y in range(0, size):    # iterate through b's rows
            for z in range(0, size):#  iterate through y's row
                return_matrix[x][y] += a[x][z] * b[z][y]
    return return_matrix

N = 4000
# Set up total number of particles that need to generate


def generate_point(n: int):
    """This function take in total number of particles and generate random
    height, z, and calculate angle, theta. Then store them into two separated
    arrays and return those two arrays.
    """
    arr_z = []
    arr_theta = []
    for _ in range(n):
        z = random.uniform(-1.0, 1.0)
        arr_z.append(z)
        if z >= 0:
            angle = 90 + np.arctan(np.sqrt(1 - z ** 2) / z) * 180 / np.pi
        else:
            angle = 270 + np.arctan(np.sqrt(1 - z ** 2) / z) * 180 / np.pi
        arr_theta.append(angle)
    return arr_z, arr_theta


def relative_prob(a: Tuple[int, int], b: Tuple[int, int], arr: List):
    """This function take in two range in the form of tuple and an array of
    number. Return the relative probability of getting this two range of angle
    """
    arr_a = []
    arr_b = []
    for angle in arr:
        if a[0] < angle < a[1]:
            arr_a.append(theta)
        if b[0] < angle < b[1]:
            arr_b.append(theta)
    return len(arr_a) / len(arr_b)


z, theta = generate_point(N)
# Create the array for z and theta
plt.plot(z, theta, "bo")
plt.title('Input z verse angle of theta')
plt.xlabel('Value of z')
plt.ylabel('Angle of theta in degree')
plt.legend()
# Plot it
plt.show()
# Show the graph

print(relative_prob((170, 190), (90, 110), theta))
# Calculate the relative probability


if __name__ == "__main__":
    # solution_to_1c(0.47, 0.0, 0.0, 8.17, 0.0001, 1)
    # solution_to_1d(0.47, 0.0, 0.0, 8.17, 0.0001, 1)
    solution_to_problem_3_N()
    #solution_to_problem_3_N_cube()
