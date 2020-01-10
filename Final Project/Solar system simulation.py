import numpy as np
from matplotlib import pyplot as plt
import time

G = 6.67e-11

#AU = 1.496e11

"""
Caculates the acceleration due to gravity
"""
def gravitational_acceleration(cur_position, position, m2):
    r = np.sqrt(np.sum(np.square(cur_position - position)))
    acceleration = -G * m2 * (cur_position - position) / r ** 3
    return acceleration.transpose()

"""
Iterate through each planet and calculate the combine acceleration of that planet due to all the other planet 
"""
def Total_acc(planet, position, masses, i, total_acceleration, radius):
    data = position[planet][:, i]
    for x in position.keys():
        if x != planet:
            total_acceleration += gravitational_acceleration(data, position[x][:, i], masses[x])
            if planet == 'earth':   # This is the part where the accelration cause by the engines are account for
                r = np.sqrt(np.sum(np.square(position['earth'][:, i])))
                if radius[planet] + radius[x] >= r:     # Check if the Earth hits other planets
                    print('Collision with ', x, ' simulation failed')
                total_acceleration[0] -= 2e-5 * position['earth'][:, i][1] / r
                total_acceleration[1] -= 2e-5 * position['earth'][:, i][0] / r
    return total_acceleration

"""
Main program
"""
def simulation(simulation_length):
    # Set up all the orbital parameters of the planets
    sun_mass = 1988500e24
    sun_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    sun_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)

    mercury_mass = 0.33011e24
    mercury_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    mercury_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    mercury_position[0][0] = 69.82e9
    mercury_velocity[1][0] = 38.86e3

    venus_mass = 4.8675e24
    venus_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    venus_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    venus_position[0][0] = 108.94e9
    venus_velocity[1][0] = 34.79e3

    earth_mass = 5.9723e24
    earth_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    earth_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    earth_position[0][0] = 152.10e9
    earth_velocity[1][0] = 29.29e3

    mars_mass = 0.64171e24
    mars_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    mars_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    mars_position[0][0] = 249.23e9
    mars_velocity[1][0] = 21.97e3

    jupiter_mass = 1898.19e24
    jupiter_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    jupiter_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    jupiter_position[0][0] = 816.62e9
    jupiter_velocity[1][0] = 12.44e3

    saturn_mass = 568.34e24
    saturn_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    saturn_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    saturn_position[0][0] = 1514.50e9
    saturn_velocity[1][0] = 9.09e3

    uranus_mass = 86.813e24
    uranus_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    uranus_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    uranus_position[0][0] = 3003.62e9
    uranus_velocity[1][0] = 6.49e3

    neptune_mass = 102.413e24
    neptune_position = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    neptune_velocity = np.zeros(simulation_length * 2 + 2, dtype=np.longdouble).reshape(2, simulation_length + 1)
    neptune_position[0][0] = 4545.67e9
    neptune_velocity[1][0] = 5.37e3

    # Put all the information of the planet into a dictionary for easy access
    position = {'mercury': mercury_position, 'venus': venus_position, 'earth': earth_position,
                'mars': mars_position, 'jupiter': jupiter_position, 'saturn': saturn_position,
                'uranus': uranus_position, 'neptune': neptune_position, 'sun': sun_position}

    mass = {'mercury': mercury_mass, 'venus': venus_mass, 'earth': earth_mass, 'mars': mars_mass,
            'jupiter': jupiter_mass, 'saturn': saturn_mass, 'uranus': uranus_mass, 'neptune': neptune_mass,
            'sun': sun_mass}

    radius = {'mercury': 2439.7e3, 'venus': 6051.8e3, 'earth': 6378.137e3, 'mars': 3396.2e3,
            'jupiter': 69991e3, 'saturn': 60268e3, 'uranus': 25559e3, 'neptune': 24764e3,
            'sun': 695510e3}

    velocity = {'mercury': mercury_velocity, 'venus': venus_velocity, 'earth': earth_velocity,
                'mars': mars_velocity, 'jupiter': jupiter_velocity, 'saturn': saturn_velocity,
                'uranus': uranus_velocity, 'neptune': neptune_velocity, 'sun': sun_velocity}

    energy = np.zeros(simulation_length, dtype=float)
    a = time.time()  # Time the performance of the code
    # The loop below is a velet method
    for i in range(0, simulation_length):
        totalE = 0
        halfacc = {}
        # The loop below update the posistion vector of all the planets
        for planet in mass.keys():
            if planet != "sun":
                if planet == 'earth':
                    totalE += 0.5 * mass[planet] * np.sum(np.square(velocity[planet][:, i]))

                net_acc = np.zeros(2, dtype=float)
                acc = Total_acc(planet, position, mass, i, net_acc, radius)
                halfacc[planet] = acc
                position[planet][:, i + 1] = position[planet][:, i] + 3600 * (velocity[planet][:, i] + 0.5 * 3600 * acc)
        # This loop update the velocity vectors for all the planet
        for planets in mass.keys():
            if planets != "sun":
                acceleration = np.zeros(2, dtype=float)
                totalacc = Total_acc(planets, position, mass, i + 1, acceleration, radius)
                velocity[planets][:, i + 1] = velocity[planets][:, i] + 0.5 * 3600 * halfacc[planets] + 0.5 * (3600 * totalacc)

        energy[i] = totalE
    c = time.time()
    print(c - a)
    # plot all the trajectory
    for p in mass.keys():
        if p == 'sun':
            plt.plot(position[p][0], position[p][1], 'ro', label=p)
        else:
            plt.plot(position[p][0], position[p][1], label=p)
    plt.legend()
    plt.title("Earth trajectory after acceleration a = 2e-5 for 17 years")
    plt.xlabel("distance, unit: m")
    plt.ylabel("distance, unit: m")
    plt.savefig("wandering_earth_2e-5.png")
    plt.show()
    plt.figure()
    plt.plot(np.arange(0, simulation_length)/(365*24), energy, label='Earth Kinetic Energy')
    plt.plot(np.arange(0, simulation_length) / (365 * 24),
             np.full(simulation_length, 0.5 * mass['earth'] * 42.1e3 ** 2), label= 'Escape Energy')
    plt.legend()
    plt.xlabel("Years, unit: year")
    plt.ylabel("Total KE, unit: J")
    plt.title("Total KE of the earth")
    plt.savefig("wandering_earth_KE_real_2e-5.png")
    plt.show()


if __name__ == "__main__":
    simulation(365*24*17)
