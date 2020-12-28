'''
ion_simulation.py
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import random

#from rims import *
'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''
diffusion_coefficient_dict = {
    # Deff dictionary. [m^2/sec]
    "Lead Pb+2"          : 0.945,
    "Potassium K+"     : 1.960,
    "Calcium Ca+2"       : 0.793,
    "Sodium Na+"        : 1.330
}

typical_ion_speed = {
    # speed dictionary. [m^2/(V*sec)]
    "Lead Pb+2": 5.19 * pow(10,-8),
    "Potassium K+": 7.62 * pow(10,-8),
    "Calcium Ca+2": 5.19 * pow(10,-8),
    "Sodium Na+": 5.19 * pow(10,-8)
}

INTERVALS_FLASH_RATIO = 10
BOLTZMANN_CONSTANT = 1.38064852 * pow(10,-23)
AVOGADRO_NUMBER = 6.0221409 * pow(np.e,23)
TEMPERATURE = 298

'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''

def get_potential_profile(a1, a2, L, x):
    # format of elecrtron experiment, units of Vq
    return a1*np.sin(2*np.pi*x / L) + a2*np.sin(4*np.pi*x / L)

def create_arena_E():
    x=np.linspace(0,0.8) # micrometer
    return plt.plot(x,get_potential_profile(-0.25, -0.05, 0.8, x))

class ion:
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode):

        self.ion = ion
        self.diffusion = diffusion_coefficient_dict[ion]
        self.potential_profile_list= potential_profile
        self.flash_frequency = flash_frequency
        self.flash_period = 1 / flash_frequency
        self.flash_mode = flash_mode

        self.L = potential_profile[0] + potential_profile[1]
        self.loc = random.uniform(0,self.L)
        self.speed = random.uniform(-typical_ion_speed[ion],typical_ion_speed[ion])

        self.intervals_count = 0
        self.points = 100
        self.arena_count = 0
        self.loc_hist = 0

    def get_intervals(self):
        self.interval = self.flash_period / INTERVALS_FLASH_RATIO

    def get_gamma(self):
        self.gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion

    def electric_force(self, x):
        electric_field = derivative(self.arena * ion.ratchet_mode(self), x, dx=1e-3)
        return - self.gamma * electric_field

    def noise(self):
        var = np.divide(1,np.sqrt(4 * np.pi * self.diffusion * self.interval))
        ksai = np.random.normal(0, np.power(var,2))
        return np.multiply(ksai,np.sqrt(2 * self.diffusion * self.interval))

    def create_arena(self):
        a = self.potential_profile_list[0]
        b = self.potential_profile_list[1]
        A = self.potential_profile_list[2]
        x = np.linspace(0, a+b)
        f1=A * np.divide(x,a)
        f2=A * np.divide((x-(a+b)),(-b))
        step = 0.5*(np.sign(x-a) +1)
        f=  f1 -step*f1 + step* f2
        self.arena = f
        # plt.plot(x, f)
        # plt.xlabel("X [um]")
        # plt.ylabel("V [v]")

    def ratchet_mode(self):
        mode = np.divide(self.intervals_count ,self.flash_period)
        if mode % 2:
            return self.flash_mode
        else:
            return 1

    def get_new_x(self):
        xt = self.loc
        new_x = xt + np.multiply(self.electric_force(xt), self.interval) + self.noise

        while new_x > self.L:
            new_x -= self.L
            self.arena_count += 1

        while new_x < self.L:
            new_x += self.L
            self.arena_count -= 1

        return new_x

    def simulate_ion(self):
        ion.create_arena(self)
        new_x =0
        while self.intervals_count <= self.points:

            new_x = ion.get_new_x(self)
            self.loc = new_x
            self.intervals_count += 1

        return self.L * self.arena_count + new_x


    #def get_max_x0(self, ):
    


