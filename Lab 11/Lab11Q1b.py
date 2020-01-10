import numpy as np
from random import random, seed
from matplotlib import pyplot as plt
tau = 1e4



# Generate a gaussian distribution
def guassian():
    r = np.sqrt(-2 * np.log(1 - random()))
    theta = 2 * np.pi * random()
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x, y

# This Function return the result for eq18
def Eq_18(x, y):
    return np.cos(x) + np.cos(np.sqrt(2) * x) + np.cos(np.sqrt(3) * x) + (y - 1)**2


def Question_1b(T_max, T_min):
    t = 0
    T = T_max
    x_array = []
    y_array = []
    x = 2
    y = 2
    old_value = Eq_18(x, y)
    while T > T_min:
        t += 1
        T = T_max * np.exp(-t / tau)

        dx, dy = guassian()
        while not (0 < x + dx < 50 and -20 < y + dy < 20):
            # If x + dx or y + dy is not in the desognated range, recalculate dx and dy.
            dx, dy = guassian()
        new_x = x + dx
        new_y = y + dy

        new_value = Eq_18(new_x, new_y)
        delta = new_value - old_value

        if random() <= np.exp(-delta / T):
            x = new_x
            y = new_y
            old_value = new_value
            x_array.append(new_x)
            y_array.append(new_y)
    print(x, y)
    print(t)
    plt.figure()
    plt.plot(np.asarray(x_array) - 2, "bo")
    plt.xlabel("Simulation iteration")
    plt.ylabel("Distance")
    plt.title("x distance between the new point and x = 2")
    plt.savefig("X_distance")
    plt.figure()
    plt.plot(np.asarray(y_array) - 1, "bo")
    plt.xlabel("Simulation iteration")
    plt.ylabel("Distance")
    plt.title("x distance between the new point and x = 2")
    plt.savefig("Y_distance")

if __name__ =="__main__":
    Question_1b(10000, 1e-7)
