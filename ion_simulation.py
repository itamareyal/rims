'''
ion_simulation.py
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''
import sympy as sym
from sympy import *
from sympy import Heaviside, S
import numpy as np
import random
import matplotlib.pyplot as plt

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
BOLTZMANN_CONSTANT = 1.38064852 * pow(10,-23)
AVOGADRO_NUMBER = 6.0221409 * pow(np.e,23)
TEMPERATURE = 298

'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''

class ion:
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode):

        self.ion = ion
        self.diffusion = diffusion_coefficient_dict[ion]
        self.potential_profile_list= potential_profile
        self.flash_frequency = flash_frequency
        self.flash_period = 1 / flash_frequency
        self.flash_mode = flash_mode

        if potential_profile[3] == 1:
            self.L = potential_profile[0]
        else:
            self.L = potential_profile[0] + potential_profile[1]
        self.loc = random.uniform(0,self.L)

        self.intervals_count = 0
        self.points = 100
        self.arena_count = 0
        self.loc_hist = 0

        # attributes calculated in simulation
        self.arena = -1
        self.potential_profile = -1
        self.interval = -1
        self.gamma = -1


    def get_intervals(self):
        self.interval = self.flash_period / INTERVALS_FLASH_RATIO


    def get_gamma(self):
        self.gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion


    def electric_force(self, x):
        electric_field = self.electric_field[x]
        return  self.gamma * electric_field


    def noise(self):
        #var = np.divide(1,np.sqrt(4 * np.pi * self.diffusion * self.interval))
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




    def electric_field(self):
        # Description: derives the electric field from the potential, E(x).
        # Parameters: self
        # Return: saves E(x) as attribute self.electric field and V(x) as self.potential_profile.
        #x = Symbol('x')

        if self.potential_profile_list[3] == 2: #sin
            L = self.potential_profile_list[0]
            a1 = self.potential_profile_list[1]
            a2 = self.potential_profile_list[2]
            # teta1 = np.divide(np.multiply(x , np.multiply(2,np.pi)),L)
            # teta2 = np.divide(np.multiply(x, np.multiply(4, np.pi)), L)
            #V= a1 * sym.sin(2*sym.pi *x / L) + a2 * sym.sin(4*sym.pi *x / L)
            #E = V.diff(x)
            x = np.linspace(0, L)
            V = a1 * np.sin(2*np.pi *x / L) + a2 * np.sin(4*np.pi *x / L)
            E = -np.gradient(V)
            # plt.plot(V)
            # plt.plot(E)
            # plt.show()

        else: #saw
            a = self.potential_profile_list[0]
            b = self.potential_profile_list[1]
            A = self.potential_profile_list[2]

            x = np.linspace(0,a+b)
            f1=A * np.divide(x,a)
            #e1 = f1.diff(x)
            f2=A * np.divide((x-(a+b)),(-b))
            #e2 = f2.diff(x)
            #step = 0.5*(np.sign(x-a) +1)
            step = np.heaviside(x-a,1)
            V=  f1 -step*f1 + step* f2
            E= -np.gradient(V)
            # plt.plot(V)
            # plt.plot(E)
            # plt.show()


        self.electric_field = E
        self.potential_profile = V
        return


    def ratchet_mode(self):
        # Description: checks if the ratchet is on or off (or negative).
        # Parameters: self
        # Return: 1 for on, 0 for off, -1 for negative.
        mode = np.divide(self.intervals_count ,self.flash_period)
        if mode % 2:
            return self.flash_mode
        else:
            return 1


    def get_new_x(self):
        # Description: calculate location for next iteration, update to arena count.
        # Parameters: self
        # Return: new location of the ion.
        xt = self.loc
        new_x = xt + np.multiply(self.electric_force(xt), self.interval) + self.noise()

        while new_x > self.L:
            new_x -= self.L
            self.arena_count += 1
        while new_x < self.L:
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
        ion.electric_field(self)
        new_x =0
        while self.intervals_count <= self.points:

            new_x = ion.get_new_x(self)
            self.loc = new_x
            self.intervals_count += 1

        return self.L * self.arena_count + new_x

