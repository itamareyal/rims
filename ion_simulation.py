'''
ion_simulation.py

Calculates the location of an ion for every time interval.
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import random
from defines import *


'''----------------------------------------------------------------------
                            IMPLEMENTATIONS
----------------------------------------------------------------------'''

class ion:
    def __init__(self, ion_subject, L, t_vec, E, path):
        """
        Instance holding ion under test
        """
        '''ratchet attributes'''
        self.diffusion = diffusion_coefficient_dict[ion_subject]
        self.electric_field_mat = E
        self.flash_frequency = 1 / t_vec[-1]
        self.flash_period = t_vec[-1]
        self.time_vec = t_vec
        self.L = L

        '''simulation attributes'''
        self.ion = ion_subject
        self.loc = random.uniform(0, self.L)
        self.x0 = self.loc
        self.intervals_count = 0
        self.points = POINTS
        self.arena_count = 0
        self.path = path

        '''result attributes'''
        self.velocity = 0
        self.arena = -1
        self.interval = -1
        self.gamma = -1
        self.vi = -1

    def get_intervals(self):
        self.interval = ((1 / INTERVALS_FLASH_RATIO) * self.L) ** 2 / (2 * self.diffusion)

    def get_gamma(self):
        self.gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion

    def electric_force(self, x):
        index_in_array = int(x * self.electric_field_mat.shape[1] / self.L)
        mode = self.ratchet_mode()
        electric_field = self.electric_field_mat[mode][index_in_array]
        return self.gamma * electric_field

    def noise(self):
        ksai = np.random.normal(0, 1)
        return np.multiply(ksai, np.sqrt(2 * self.diffusion * self.interval))

    def ratchet_mode(self):
        """
        :return: the index of the active potential profile based on the current time
        """
        t_prime = np.mod(self.intervals_count * self.interval, self.flash_period)
        mode = 0
        while t_prime > self.time_vec[mode]:
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
        ion.get_intervals(self)
        ion.get_gamma(self)

        new_x = 0
        dummy_ss = int(0.97 * self.points)
        vi_list = []

        while self.intervals_count <= self.points:

            prev_loc = self.arena_count * self.L + self.loc
            new_x = ion.get_new_x(self)
            self.loc = new_x
            if self.intervals_count > dummy_ss:
                vi = (new_x * self.arena_count * self.L - prev_loc) / self.interval
                vi_list.append(vi)
            self.intervals_count += 1

        vi_vector = np.array(vi_list)
        self.vi = np.average(vi_vector)
        ret_val = self.L * self.arena_count + new_x

        '''keeping final location in range [0, num of ratchets times arena size] to get steady state'''
        while ret_val > RATCHETS_IN_SYSTEM * self.L:
            ret_val -= RATCHETS_IN_SYSTEM * self.L
        while ret_val < 0:
            ret_val += RATCHETS_IN_SYSTEM * self.L

        return ret_val

