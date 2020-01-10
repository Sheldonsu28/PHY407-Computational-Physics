import numpy as np
import math as m
import pylab as plt
from scipy import constants as con
from guassxw import gaussxw

h = 6.62607004 * 10 ** -34
c = con.speed_of_light
kb = con.Boltzmann


def I_Lamb_2a(L, t):
    numerator = L ** -5
    denominator = (m.e ** ((h * c) / (L * kb * t))) - 1
    return numerator / denominator


def I_x_2a(x, t):
    numerator = (x / (x + 1)) ** -5
    denominator = ((m.e ** ((h * c * (x + 1)) / (x * kb * t))) - 1) * (x + 1) ** 2
    return numerator / denominator


def Q2_Gaussian_integration_efficency(temperature, x, w, x2, w2):
    I_x = []
    I_lambda = []

    lambda_x, lambda_w = x, w
    lambda_x = np.asarray(lambda_x)
    lambda_w = np.asarray(lambda_w)

    lambda_x_partition = 0.5 * (400 * 10 ** -9) * lambda_x + 0.5 * (380 + 780) * 10 ** -9

    lambda_w_partition = 0.5 * (400 * 10 ** -9) * lambda_w

    x_x, x_w = x2, w2
    x_x_partition = 0.5 * 0.0001 * x_x + 0.5 * 0.0001

    x_w_partition = 0.5 * 0.0001 * x_w

    for x in range(0, 4000):
        I_x.append(x_w_partition[x] * I_x_2a(x_x_partition[x], temperature))

    for y in range(0, 400):
        I_lambda.append(lambda_w_partition[y] * I_Lamb_2a(lambda_x_partition[y], temperature))

    return sum(I_lambda) / sum(I_x)


def Question_2a():
    temperature = np.arange(300, 10000, 1)
    efficency = np.zeros(10000 - 300, dtype=float)
    count = 0
    x, w = gaussxw(400)
    x2, w2 = gaussxw(4000)
    for t in temperature:
        efficency[count] = Q2_Gaussian_integration_efficency(t, x, w, x2, w2)
        count += 1
    plt.plot(temperature, efficency, "bo")

    plt.ylabel("Efficency of the light bulb")
    plt.xlabel("Temperature in K")
    plt.title("The efficency of the light bulb vs temperature")
    plt.savefig("Question 2a.png")
    plt.show()


def derivitive_Question_2_b(x, w, x2, w2, temperature, diff):
    return (Q2_Gaussian_integration_efficency(temperature + diff, x, w, x2, w2) - Q2_Gaussian_integration_efficency(
        temperature, x, w, x2, w2)) / diff


def Question2_b_helper(k):
    x, w = gaussxw(400)
    x2, w2 = gaussxw(4000)
    accuracy = 0.5
    t_1, t_2 = k, 10000
    error = 2
    while error > accuracy:
        midpoint = (t_1 + t_2) / 2
        point_1 = derivitive_Question_2_b(x, w, x2, w2, midpoint, 0.001)
        point_2 = derivitive_Question_2_b(x, w, x2, w2, t_2, 0.001)
        if point_1 > 0 and point_2 > 0:
            t_1 = midpoint
        else:
            t_2 = midpoint
        error = abs(t_1 - t_2)
    return (t_1 + t_2) / 2


def Question2_b_solution():
    inital = Question2_b_helper(6000)
    after = Question2_b_helper(inital)
    while abs(after - inital) > 1:
        after = Question2_b_helper(inital)
        inital = Question2_b_helper(after)
    print(inital)


def Question3_a():
    y = []
    x = []
    threshhold = 1.0e-6
    C = 0
    for i in range(0, 300):
        x.append(C)
        initial_value = 1
        error = 1.0
        new_value = 0
        while error > threshhold:  # Question 3a
            initial_value, new_value = 1 - m.e ** (-C * initial_value), initial_value
            error = abs((new_value - initial_value) * -C * m.e ** (C * initial_value))
        y.append(initial_value)
        C += 0.01
    plt.plot(x, y)
    plt.xlabel("Value of C")
    plt.ylabel("Solution of x")
    plt.title("Values of C vs solution of x")
    plt.savefig("Question3a.png")
    plt.show()


def Question3_b():
    threshhold = 1.0e-6
    C = 2

    i_value = 1
    error = 1.0
    new_value = 0
    count_relaxation = 0
    while error > threshhold:  # Question 3 b (6.11 part b)
        i_value, new_value = 1 - m.e ** (-C * i_value), i_value
        error = abs((new_value - i_value) * -C * m.e ** (C * i_value))
        count_relaxation += 1

    w = 0.5

    initial_value = 1
    error = 1.0
    new_value = 0
    count_overrelaxation = 0
    while error > threshhold:  # Overrlaxation method
        initial_value, new_value = (1 + w) * (1 - m.e ** (-C * initial_value)) * w * initial_value, initial_value
        error = abs((new_value - initial_value) * -C * m.e ** (C * initial_value))
        count_overrelaxation += 1

    print("Iteration for Relaxation method when C = 2:", count_relaxation)
    print("Iteration for OverRelaxation method when C = 2:", count_overrelaxation)


def Question3_c():
    threshhold = 1.0e-6
    x_inital = 2
    y_inital = 2
    x_new = 0
    y_new = 0
    error_x = 1.0
    error_y = 1.0
    count_x = 0
    count_y = 0
    a = 1
    b = 2
    while error_x > threshhold and count_x < 5:
        x_inital, x_new = y_inital * (a + x_inital ** 2), x_inital
        error_x = abs((x_new - x_inital) * y_inital * (a + 2 * x_inital))
        count_x += 1

    while error_y > threshhold and count_y < 5:
        y_inital, y_new = b / (a + x_inital ** 2), y_inital
        error_y = abs((y_new - y_inital) * 2 * b * x_inital / ((a + x_inital ** 2) ** 2))
        count_y += 1

    if count_x >= 5 or count_y >= 5:
        print("Equations failed to converage, solution not found")
    else:
        print("x solution: ", x_inital)
        print("y solution: ", y_inital)


def Question3_c_part_c():
    threshhold = 1.0e-6
    x_inital = 2
    y_inital = 0.4
    x_new = 0
    y_new = 0
    error_x = 1.0
    error_y = 1.0
    count_x = 0
    count_y = 0
    a = 1
    b = 2
    while error_x > threshhold and count_x < 1000000:
        x_inital, x_new = m.sqrt((x_inital / y_inital) + a), x_inital
        error_x = abs((x_new - x_inital) * (2 * y_inital * m.sqrt((x_inital / y_inital) + a)) ** -1)
        count_x += 1

    while error_y > threshhold and count_y < 1000000:
        y_inital, y_new = (b - a * y_inital) / x_inital ** 2, y_inital
        error_y = abs((y_new - y_inital) * -a / (x_inital ** 2))
        count_y += 1

    if count_x >= 1000000 and count_y >= 1000000:
        print("Equations failed to converage, solution not found")
    else:
        print("x solution: ", x_inital)
        print("y solution: ", y_inital)


if __name__ == "__main__":
    # Question_2a()
    # Question2_b_solution()
    # Question3_a()
    # Question3_b()
    Question3_c_part_c()
