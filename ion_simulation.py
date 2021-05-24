'''
ion_simulation.py

Calculates the location of an ion for every time interval.
'''


'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import random
from defines import *
from current_calc import *

'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''

class ion:
    def __init__(self, rims):
        """
        Instance holding ion under test
        """
        '''ratchet attributes'''
        self.diffusion = rims.diffusion
        self.electric_field_mat = rims.electric_field_mat
        self.time_vec = rims.time_vec
        self.flash_frequency = rims.flash_frequency
        self.flash_period = rims.flash_period
        self.intervals_in_period = rims.intervals_in_period
        self.L = rims.L

        '''simulation attributes'''
        self.ion = rims.ion
        self.loc = random.uniform(0, self.L)
        self.x0 = self.loc
        self.intervals_count = 0
        self.points = POINTS
        self.arena_count = 0
        self.path = rims.path_for_output

        '''result attributes'''
        self.absolute_final_loc = None
        self.velocity = 0
        self.interval = rims.interval
        self.gamma = rims.gamma

    def electric_force(self, x):
        index_in_array = int(x * self.electric_field_mat.shape[1] / self.L)
        mode = self.ratchet_mode()
        electric_field = self.electric_field_mat[mode][index_in_array]
        return electric_field / self.gamma

    def noise(self):
        ksai = np.random.normal(0, 1)
        return np.multiply(ksai, np.sqrt(2 * self.diffusion * self.interval))

    def ratchet_mode(self):
        """
        :return: the index of the active potential profile based on the current time
        """
        t_prime = np.mod(self.intervals_count * self.interval, self.flash_period)
        mode = 0
        while t_prime >= self.time_vec[mode]:
            mode += 1
        return mode

    def get_new_x(self):
        """
        calculate location for next iteration, update to arena count.
        :return: new location of the ion.
        """
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
        """
        simulate ion movement over number of iterations.
        :return: location of the ion at the nd of the simulation.
        """
        self.intervals_count = 0
        relative_x0 = (self.arena_count * self.L) + self.loc
        while self.intervals_count < self.intervals_in_period:
            # prev_loc = (self.arena_count * self.L) + self.loc
            self.loc = ion.get_new_x(self)
            # next_loc = (self.arena_count * self.L) + self.loc
            # self.velocity_list.append(calculate_v(next_loc, prev_loc, self.interval))           # Each ion creates a velocity_list which includes all of the velocities sampled during the simulation
            self.intervals_count += 1

        self.absolute_final_loc = self.L * self.arena_count + self.loc
        self.velocity = (self.absolute_final_loc - relative_x0) * self.flash_frequency
        return




