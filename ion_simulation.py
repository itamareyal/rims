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
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode, dc, E, V, path):

        self.ion = ion
        self.diffusion = diffusion_coefficient_dict[ion]
        self.potential_profile_list= potential_profile

        self.electric_field = E
        self.potential_profile = V
        self.velocity = 0

        self.flash_frequency = flash_frequency
        self.flash_period = 1 / flash_frequency
        self.flash_mode = flash_mode
        self.dc = dc

        if potential_profile[3] == 2:       # sin
            self.L = potential_profile[0]
        else:
            self.L = potential_profile[0] + potential_profile[1]
        self.loc = random.uniform(0, self.L)
        self.x0= self.loc

        self.intervals_count = 0
        self.points = POINTS
        self.arena_count = 0
        self.path = path

        # attributes calculated in simulation
        self.arena = -1
        self.potential_profile = -1
        self.interval = -1
        self.num_of_intervals_for_current = -1
        self.gamma = -1
        self.steady_state_velocity = -1
        self.steady_state = False
        self.velocity_list = []

    def get_num_of_intervals_for_current(self):
        self.num_of_intervals_for_current = int(self.flash_period / self.interval)

    def get_intervals(self):
        self.interval = ((1 / INTERVALS_FLASH_RATIO) * self.L) ** 2 / (2 * self.diffusion)

    def get_gamma(self):
        self.gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion

    def electric_force(self, x):
        index_in_array = int(x * RESOLUTION / self.L)
        electric_field = self.electric_field[index_in_array]
        return self.gamma * electric_field * self.ratchet_mode()

    def noise(self):
        ksai = np.random.normal(0, 1)
        return np.multiply(ksai, np.sqrt(2 * self.diffusion * self.interval))

    def create_arena(self):
        """
        Creates the potential profile V(x) over one period. for plotting.
        saves the arena as attribute self.arena
        """
        if self.potential_profile_list[3] == 2:         # sin
            L = self.potential_profile_list[0]
            a1 = self.potential_profile_list[1]
            a2 = self.potential_profile_list[2]
            x = np.linspace(0, L)
            self.arena = a1 * np.sin(2 * np.pi * x / L) + a2 * np.sin(4 * np.pi * x / L)

        else:                                           # saw
            a = self.potential_profile_list[0]
            b = self.potential_profile_list[1]
            A = self.potential_profile_list[2]
            x = np.linspace(0, a+b)
            f1 = A * np.divide(x, a)
            f2 = A * np.divide((x-(a+b)), (-b))
            step = 0.5*(np.sign(x-a) + 1)
            f = f1 - step * f1 + step * f2
            self.arena = f

    def ratchet_mode(self):
        """
        checks if the ratchet is on or off (or negative).
        :return: 1 for on, 0 for off, -1 for negative.
        """
        t_prime = np.mod(self.intervals_count * self.interval, self.flash_period)
        if t_prime < self.dc * self.flash_period:
            return 1
        else:
            return self.flash_mode

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

    def calculate_velocity(self, v_list):
        n = self.num_of_intervals_for_current
        v_array = v_list[-n:]
        np.array(v_array)
        v = np.average(v_array)
        self.velocity_list.append(v)

    # def check_for_steady_state(self, vt_list):
    #     if len(vt_list) < 5:
    #         return
    #     last_vi_margin = vt_list[-1] / 10
    #     for i in range(2, 6, 1):
    #         if (vt_list[-i] <= vt_list[-1] - last_vi_margin) or (vt_list[-i] >= vt_list[-1] + last_vi_margin):
    #             return
    #     self.steady_state = True
    #     return

    def simulate_ion(self, number_of_cycles_per_ion):
        """
        simulate ion movement over number of iterations.
        :return: location of the ion at the nd of the simulation.
        """
        ion.get_intervals(self)
        ion.get_num_of_intervals_for_current(self)
        ion.get_gamma(self)
        ion.create_arena(self)

        vi_list = []

        while self.intervals_count <= (number_of_cycles_per_ion * self.num_of_intervals_for_current) and self.intervals_count <= self.points:
            old_location = (self.arena_count * self.L) + self.loc
            self.loc = ion.get_new_x(self)
            new_location = (self.arena_count * self.L) + self.loc
            self.velocity_list.append(calculate_v(new_location, old_location, self.interval))           # Each ion creates a velocity_list which includes all of the velocities sampled during the simulation
            self.intervals_count += 1

            # if (self.intervals_count % self.num_of_intervals_for_current) == 0:
            #     ion.calculate_velocity(self, vi_list)

            # ion.check_for_steady_state(self, self.velocity_list)
            # if self.steady_state:
            #     break

        # print("Reached Steady State after " + str(self.intervals_count + 1) + " intervals")
        # vt_vector = np.array(self.velocity_list[-5:])
        # self.steady_state_velocity = np.average(vt_vector)

        ret_val = self.L * self.arena_count + self.loc

        # keeping final location in range [0, num of ratchets times arena size] to get steady state
        while ret_val > RATCHETS_IN_SYSTEM * self.L:
            ret_val -= RATCHETS_IN_SYSTEM * self.L
        while ret_val < 0:
            ret_val += RATCHETS_IN_SYSTEM * self.L
        # calculate group number
        return ret_val



