import numpy as np
import pylab as plt
import time


def Question1(omega):
    height = 80
    width = 200
    phi = np.zeros([height + 2, width + 2], float)
    shape = np.ones([height + 2, width + 2], int)
    AH = 20/80
    AB = 5/90
    CB = 20/40
    for a in range(92, 111, 1):             # initalize all the edges
        for b in range(1, 41, 1):
            shape[b][a] = 0

    for width in range(0, 202, 1):
        shape[0][width] = 0
        shape[81][width] = 0
        phi[80][width] = 20
        shape[80][width] = 0

    for height in range(0, 82, 1):
        shape[height][0] = 0
        shape[height][201] = 0
        phi[height][1] = AH*(height - 1)
        phi[height][200] = AH * (height - 1)
        shape[height][1] = 0
        shape[height][200] = 0

    for x in range(1, 91, 1):
        phi[1][x] = AB * (x - 1)
        phi[1][x + 111] = -AB * (x - 1) + 5
        shape[1][x] = 0
        shape[1][x + 111] = 0
    for x in range(91, 111):
        phi[40][x] = 25
        shape[40][x] = 0
    for y in range(1, 41):
        phi[y][91] = CB * y
        shape[y][91] = 0
        phi[y][111] = CB * y
        shape[y][111] = 0



    phi[1][1] = 0
    phi[1][91] = 5
    phi[41][91] = 25
    phi[42][111] = 25
    phi[1][111] = 5
    phi[1][200] = 0
    phi[80][200] = 20
    phi[80][1] = 20

    shape[1][1] = 1
    shape[1][91] = 1
    shape[41][91] = 1
    shape[42][111] = 1
    shape[1][111] = 1
    shape[1][200] = 1
    shape[80][200] = 1
    shape[80][1] = 1


    error = 1
    while error > 1e-6:                         # main loop
        old = np.copy(phi)
        for x in range(0, width + 1):
            for y in range(0, height + 1):
                if shape[y][x] == 1:
                    phi[y][x] = ((1 + omega) / 4) * (
                            phi[y + 1, x] + phi[y - 1, x] + phi[y, x + 1] + phi[y, x - 1]) - omega * phi[y][x]
        error = np.max(np.abs(phi - old))
    print(phi[10][25])
    plt.imshow(phi, cmap=plt.set_cmap("coolwarm"), extent=(0, 20, 8, 0))
    plt.colorbar(shrink=0.6, label="Temperature")
    plt.xlabel("Unit: cm")
    plt.ylabel("Unit cm")
    plt.title("Temperature with Omega = 0.9")

    plt.savefig("Quetion 1 9e-1.png")
    plt.show()


def Helper_Q2(t):
    return 10 + 12 * np.sin((2 * np.pi * t) / 365)


def Question2():
    height = 200
    spacing = 0.1
    step = 0.01
    results = []
    D = 0.1
    c = step * D / (spacing * spacing)
    ground = np.full(height, 10, float)
    ground[height - 1] = 11
    ground[0] = Helper_Q2(0)
    results.append(ground)

    for t in range(0, 3650, 1):
        ground = np.copy(results[-1])
        ground[0] = Helper_Q2(t)
        for i in range(0, 100):
            for point in range(1, height - 1):
                ground[point] = ground[point] + c * (ground[point + 1] + ground[point - 1] - 2 * ground[point])
        results.append(ground)

    depth = np.arange(0, 20, 0.1)

    plt.plot(depth, results[3285], label="First three months")
    plt.plot(depth, results[3285 + 93], label="Second three months")
    plt.plot(depth, results[3285 + 2 * 93], label="Third three months")
    plt.plot(depth, results[3580], label="Last three months")
    plt.legend()
    plt.ylim((-3, 23))
    plt.xlabel("Depth unit: m")
    plt.ylabel("Temperature unit: Celcius")
    plt.title("Temperture vs depth")
    plt.savefig("Q2.png")
    plt.show()










def Question3(s):
    epsilon = 1
    dx = 0.02
    dt = 0.005
    Lx = 2*np.pi
    Tf = 2
    beta = epsilon * dt/dx

    x = np.arange(0, Lx, dx)
    time_1 = np.arange(0, Tf, dt)

    u = np.zeros([len(x), len(time_1)])

    for x_value in range(1, len(x) - 1):                            # Forward method
        u[x_value, 0] = np.sin(x[x_value])
        u[x_value, 1] = u[x_value, 0] + dt * np.cos(x[x_value])

    for t in range(1, len(time_1) - 1):
        for x_1 in range(1, len(x) - 1):
            u[x_1, t + 1] = u[x_1, t - 1] - beta / 2 * (u[x_1 + 1, t] ** 2 - u[x_1 - 1, t] ** 2)
    plt.plot(u[:,int(s * 200)])
    plt.title('Wave Function at time = 1.5')
    plt.xlabel('X Position')
    plt.ylabel('Wave Velocity')
    plt.savefig("Q30t1e5.png")
    plt.show()


if __name__ =="__main__":
    #Question1(0.9)
    #Question2()
    Question3(1.5)
