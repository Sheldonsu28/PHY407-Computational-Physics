import numpy as np
from matplotlib import pyplot as plt
import time
from multiprocessing import Process, Queue
import multiprocessing as mp
import ctypes

G = 6.67e-11


def gravitational_acceleration(cur_position, position, m2):
    r = np.sqrt(np.sum(np.square(cur_position - position)))
    acceleration = -G * m2 * (cur_position - position) / r ** 3
    return acceleration.transpose()


def update_posistion(p0, p1, p2, p3, p4, p5, p6, p7, p8, v0, v1, v2, v3, v4, v5, v6, v7, v8, mass, conn, planet,
                     sim_len):
    position = {'mercury': np.frombuffer(p1.get_obj()).reshape(2, sim_len + 1),
                'venus': np.frombuffer(p2.get_obj()).reshape(2, sim_len + 1),
                'earth': np.frombuffer(p3.get_obj()).reshape(2, sim_len + 1),
                'mars': np.frombuffer(p4.get_obj()).reshape(2, sim_len + 1),
                'jupiter': np.frombuffer(p5.get_obj()).reshape(2, sim_len + 1),
                'saturn': np.frombuffer(p6.get_obj()).reshape(2, sim_len + 1),
                'uranus': np.frombuffer(p7.get_obj()).reshape(2, sim_len + 1),
                'neptune': np.frombuffer(p8.get_obj()).reshape(2, sim_len + 1),
                'sun': np.frombuffer(p0.get_obj()).reshape(2, sim_len + 1)}

    velocity = {'mercury': np.frombuffer(v1.get_obj()).reshape(2, sim_len + 1),
                'venus': np.frombuffer(v2.get_obj()).reshape(2, sim_len + 1),
                'earth': np.frombuffer(v3.get_obj()).reshape(2, sim_len + 1),
                'mars': np.frombuffer(v4.get_obj()).reshape(2, sim_len + 1),
                'jupiter': np.frombuffer(v5.get_obj()).reshape(2, sim_len + 1),
                'saturn': np.frombuffer(v6.get_obj()).reshape(2, sim_len + 1),
                'uranus': np.frombuffer(v7.get_obj()).reshape(2, sim_len + 1),
                'neptune': np.frombuffer(v8.get_obj()).reshape(2, sim_len + 1),
                'sun': np.frombuffer(v0.get_obj()).reshape(2, sim_len + 1)}
    while True:

        i = conn.recv()

        net_acc = np.zeros(2, dtype=float)

        data = position[planet][:, i]
        for x in position.keys():
            if x != planet:
                r = np.sqrt(np.sum(np.square(data - position[x][:, i])))
                net_acc += -G * mass[x] * (data - position[x][:, i]) / r ** 3

        position[planet][:, i + 1] = velocity[planet][:, i] + 0.5 * 3608 * net_acc + 0.5 * (3600 * net_acc)
        conn.send(1)
        if i == sim_len - 1:
            break


def update_velocity(p0, p1, p2, p3, p4, p5, p6, p7, p8, v0, v1, v2, v3, v4, v5, v6, v7, v8, mass, conn, planets,
                    sim_len):
    position = {'mercury': np.frombuffer(p1.get_obj()).reshape(2, sim_len + 1),
                'venus': np.frombuffer(p2.get_obj()).reshape(2, sim_len + 1),
                'earth': np.frombuffer(p3.get_obj()).reshape(2, sim_len + 1),
                'mars': np.frombuffer(p4.get_obj()).reshape(2, sim_len + 1),
                'jupiter': np.frombuffer(p5.get_obj()).reshape(2, sim_len + 1),
                'saturn': np.frombuffer(p6.get_obj()).reshape(2, sim_len + 1),
                'uranus': np.frombuffer(p7.get_obj()).reshape(2, sim_len + 1),
                'neptune': np.frombuffer(p8.get_obj()).reshape(2, sim_len + 1),
                'sun': np.frombuffer(p0.get_obj()).reshape(2, sim_len + 1)}

    velocity = {'mercury': np.frombuffer(v1.get_obj()).reshape(2, sim_len + 1),
                'venus': np.frombuffer(v2.get_obj()).reshape(2, sim_len + 1),
                'earth': np.frombuffer(v3.get_obj()).reshape(2, sim_len + 1),
                'mars': np.frombuffer(v4.get_obj()).reshape(2, sim_len + 1),
                'jupiter': np.frombuffer(v5.get_obj()).reshape(2, sim_len + 1),
                'saturn': np.frombuffer(v6.get_obj()).reshape(2, sim_len + 1),
                'uranus': np.frombuffer(v7.get_obj()).reshape(2, sim_len + 1),
                'neptune': np.frombuffer(v8.get_obj()).reshape(2, sim_len + 1),
                'sun': np.frombuffer(v0.get_obj()).reshape(2, sim_len + 1)}

    while True:
        i = conn.recv()
        acceleration = np.zeros(2, dtype=float)

        data = position[planets][:, i]
        for x in position.keys():
            if x != planets:
                r = np.sqrt(np.sum(np.square(data - position[x][:, i])))
                acceleration += -G * mass[x] * (data - position[x][:, i]) / r ** 3

        velocity[planets][:, i + 1] = velocity[planets][:, i] + 0.5 * 3600 * acceleration + 0.5 * (3600 * acceleration)
        conn.send(1)
        if i == sim_len - 1:
            break


def simulation(simulation_length):
    sun_mass = 1988500e24
    sun_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p0 = mp.Array(ctypes.c_double, sun_position)
    sun_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v0 = mp.Array(ctypes.c_double, sun_velocity)

    mercury_mass = 0.33011e24
    mercury_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p1 = mp.Array(ctypes.c_double, mercury_position)
    mercury_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v1 = mp.Array(ctypes.c_double, mercury_velocity)
    np.frombuffer(p1.get_obj()).reshape(2, simulation_length + 1)[0][0] = 69.82e9
    np.frombuffer(v1.get_obj()).reshape(2, simulation_length + 1)[1][0] = 38.86e3

    venus_mass = 4.8675e24
    venus_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p2 = mp.Array(ctypes.c_double, venus_position)
    venus_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v2 = mp.Array(ctypes.c_double, venus_velocity)
    np.frombuffer(p2.get_obj()).reshape(2, simulation_length + 1)[0][0] = 108.94e9
    np.frombuffer(v2.get_obj()).reshape(2, simulation_length + 1)[1][0] = 34.79e3

    earth_mass = 5.9723e24
    earth_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p3 = mp.Array(ctypes.c_double, earth_position)
    earth_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v3 = mp.Array(ctypes.c_double, earth_velocity)
    np.frombuffer(p3.get_obj()).reshape(2, simulation_length + 1)[0][0] = 152.10e9
    np.frombuffer(v3.get_obj()).reshape(2, simulation_length + 1)[1][0] = 29.29e3

    mars_mass = 0.64171e24
    mars_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p4 = mp.Array(ctypes.c_double, mars_position)
    mars_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v4 = mp.Array(ctypes.c_double, mars_velocity)
    np.frombuffer(p4.get_obj()).reshape(2, simulation_length + 1)[0][0] = 249.23e9
    np.frombuffer(v4.get_obj()).reshape(2, simulation_length + 1)[1][0] = 21.97e3

    jupiter_mass = 1898.19e24
    jupiter_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p5 = mp.Array(ctypes.c_double, jupiter_position)
    jupiter_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v5 = mp.Array(ctypes.c_double, jupiter_velocity)
    np.frombuffer(p5.get_obj()).reshape(2, simulation_length + 1)[0][0] = 816.62e9
    np.frombuffer(v5.get_obj()).reshape(2, simulation_length + 1)[1][0] = 12.44e3

    saturn_mass = 568.34e24
    saturn_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p6 = mp.Array(ctypes.c_double, saturn_position)
    saturn_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v6 = mp.Array(ctypes.c_double, saturn_velocity)
    np.frombuffer(p6.get_obj()).reshape(2, simulation_length + 1)[0][0] = 1514.50e9
    np.frombuffer(v6.get_obj()).reshape(2, simulation_length + 1)[1][0] = 9.09e3

    uranus_mass = 86.813e24
    uranus_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p7 = mp.Array(ctypes.c_double, uranus_position)
    uranus_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v7 = mp.Array(ctypes.c_double, uranus_velocity)
    np.frombuffer(p7.get_obj()).reshape(2, simulation_length + 1)[0][0] = 3003.62e9
    np.frombuffer(v7.get_obj()).reshape(2, simulation_length + 1)[1][0] = 6.49e3

    neptune_mass = 102.413e24
    neptune_position = np.zeros(simulation_length * 2 + 2, dtype=float)
    p8 = mp.Array(ctypes.c_double, neptune_position)
    neptune_velocity = np.zeros(simulation_length * 2 + 2, dtype=float)
    v8 = mp.Array(ctypes.c_double, neptune_velocity)
    np.frombuffer(p8.get_obj()).reshape(2, simulation_length + 1)[0][0] = 4545.67e9
    np.frombuffer(v8.get_obj()).reshape(2, simulation_length + 1)[1][0] = 5.37e3

    position = {'mercury': np.frombuffer(p1.get_obj()).reshape(2, simulation_length + 1),
                'venus': np.frombuffer(p2.get_obj()).reshape(2, simulation_length + 1),
                'earth': np.frombuffer(p3.get_obj()).reshape(2, simulation_length + 1),
                'mars': np.frombuffer(p4.get_obj()).reshape(2, simulation_length + 1),
                'jupiter': np.frombuffer(p5.get_obj()).reshape(2, simulation_length + 1),
                'saturn': np.frombuffer(p6.get_obj()).reshape(2, simulation_length + 1),
                'uranus': np.frombuffer(p7.get_obj()).reshape(2, simulation_length + 1),
                'neptune': np.frombuffer(p8.get_obj()).reshape(2, simulation_length + 1),
                'sun': np.frombuffer(p0.get_obj()).reshape(2, simulation_length + 1)}

    mass = {'mercury': mercury_mass, 'venus': venus_mass, 'earth': earth_mass, 'mars': mars_mass,
            'jupiter': jupiter_mass, 'saturn': saturn_mass, 'uranus': uranus_mass, 'neptune': neptune_mass,
            'sun': sun_mass}

    velocity = {'mercury': mercury_velocity, 'venus': venus_velocity, 'earth': earth_velocity,
                'mars': mars_velocity, 'jupiter': jupiter_velocity, 'saturn': saturn_velocity,
                'uranus': uranus_velocity, 'neptune': neptune_velocity, 'sun': sun_velocity}

    energy = np.zeros(simulation_length + 1, dtype=float)

    position_conn_endpoint = []
    velocity_conn_endpoint = []

    for planet in mass.keys():
        if planet != "sun":
            main_p_end, child_p_end = mp.Pipe()
            position_conn_endpoint.append(main_p_end)

            main_v_end, child_v_end = mp.Pipe()
            velocity_conn_endpoint.append(main_v_end)

            processes = Process(target=update_posistion, args=(
            p0, p1, p2, p3, p4, p5, p6, p7, p8, v0, v1, v2, v3, v4, v5, v6, v7, v8, mass, child_p_end, planet,
            simulation_length))

            processes1 = Process(target=update_velocity, args=(
            p0, p1, p2, p3, p4, p5, p6, p7, p8, v0, v1, v2, v3, v4, v5, v6, v7, v8, mass, child_v_end, planet,
            simulation_length))
            processes.start()
            processes1.start()
    a = time.time()
    for i in range(1, simulation_length):
        totalE = 0
        for reciver in position_conn_endpoint:
            reciver.send(i)

        for source in position_conn_endpoint:
            source.recv()

        for reciver_v in velocity_conn_endpoint:
            reciver_v.send(i)

        for source_v in velocity_conn_endpoint:
            source_v.recv()

        energy[i] = totalE

    k = time.time()
    print(k - a)

    for p in mass.keys():
        if p == 'sun':
            plt.plot(position[p][0], position[p][1], 'ro', label=p)
        else:
            plt.plot(position[p][0], position[p][1], label=p)
    plt.legend()
    plt.show()
    plt.figure()
    plt.plot(np.arange(0, simulation_length + 1), energy)
    plt.show()


if __name__ == "__main__":
    simulation(3)
