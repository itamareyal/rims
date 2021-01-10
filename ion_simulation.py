'''
ion_simulation.py

Calculates the location of an ion for every time interval.
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import numpy as np
import random
import csv


'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''
diffusion_coefficient_dict = {
    # Deff dictionary. [m^2/sec]
    "Lead Pb+2"             : 0.945,
    "Potassium K+"          : 1.960,
    "Calcium Ca+2"          : 0.793,
    "Sodium Na+"            : 1.330,
    "Electron in Silicone"  : 0.0036
}

INTERVALS_FLASH_RATIO = 10
BOLTZMANN_CONSTANT = 8.617333262 * pow(10,-5)
AVOGADRO_NUMBER = 6.0221409 * pow(np.e,23)
TEMPERATURE = 298
RESOLUTION = 1000

'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''

class ion:
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode, dc, E, V, path):

        self.ion = ion
        self.diffusion = diffusion_coefficient_dict[ion]
        self.potential_profile_list= potential_profile

        self.electric_field = E
        self.potential_profile = V

        self.flash_frequency = flash_frequency
        self.flash_period = 1 / flash_frequency
        self.flash_mode = flash_mode
        self.dc = dc

        if potential_profile[3] == 2: #sin
            self.L = potential_profile[0]
        else:
            self.L = potential_profile[0] + potential_profile[1]
        self.loc = random.uniform(0,self.L)
        self.x0= self.loc

        self.intervals_count = 0
        self.points = 100
        self.arena_count = 0
        self.path = path

        # attributes calculated in simulation
        self.arena = -1
        self.potential_profile = -1
        self.interval = -1
        self.gamma = -1


    def get_intervals(self):
        self.interval = ((1 / INTERVALS_FLASH_RATIO) * self.L) ** 2 /(2 *self.diffusion)


    def get_gamma(self):
        self.gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion


    def electric_force(self, x):
        index_in_array = int(x * RESOLUTION / self.L)
        electric_field = self.electric_field[index_in_array]
        return  self.gamma * electric_field * self.ratchet_mode()


    def noise(self):
        ksai = np.random.normal(0, 1)
        return np.multiply(ksai,np.sqrt(2 * self.diffusion * self.interval))


    def create_arena(self):
        # Description: Creates the potential profile V(x) over one period. for plotting.
        # Parameters: self
        # Return: saves the arena as attribute self.arena

        if self.potential_profile_list[3] == 2: #sin
            L = self.potential_profile_list[0]
            a1 = self.potential_profile_list[1]
            a2 = self.potential_profile_list[2]
            x = np.linspace(0, L)
            self.arena = a1 * np.sin(2 * np.pi * x / L) + a2 * np.sin(4 * np.pi * x / L)

        else: #saw
            a = self.potential_profile_list[0]
            b = self.potential_profile_list[1]
            A = self.potential_profile_list[2]
            x = np.linspace(0, a+b)
            f1=A * np.divide(x,a)
            f2=A * np.divide((x-(a+b)),(-b))
            step = 0.5*(np.sign(x-a) +1)
            f=  f1 -step*f1 + step* f2
            self.arena = f


    def ratchet_mode(self):
        # Description: checks if the ratchet is on or off (or negative).
        # Parameters: self
        # Return: 1 for on, 0 for off, -1 for negative.

        t_prime = np.mod(self.intervals_count * self.interval ,self.flash_period)
        if t_prime < self.dc * self.flash_period:
            return 1
        else:
            return self.flash_mode


    def get_new_x(self):
        # Description: calculate location for next iteration, update to arena count.
        # Parameters: self
        # Return: new location of the ion.

        xt = self.loc
        noise = self.noise()
        fe = self.electric_force(xt)
        em = np.multiply(fe, self.interval)
        new_x = xt + em + noise

        while new_x > self.L:
            new_x -= self.L
            self.arena_count += 1
        while new_x < 0:
            new_x += self.L
            self.arena_count -= 1

        return new_x


    def simulate_ion(self):
        # Description: simulate ion movement over number of iterations.
        # Parameters: self
        # Return: location of the ion at the nd of the simulation.

        ion.get_intervals(self)
        ion.get_gamma(self)
        ion.create_arena(self)

        new_x =0
        while self.intervals_count <= self.points:

            new_x = ion.get_new_x(self)
            self.loc = new_x
            self.intervals_count += 1

        return self.L * self.arena_count + new_x

