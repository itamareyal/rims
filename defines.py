"""
defines.py

Hard coded data for the whole program
All units in cm / sec / Hz / kelvin / ev
"""

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import numpy as np


'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''

'''PARTICLE PARAMETERS'''
ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electrons in Silicon"]

diffusion_coefficient_dict = {                      # units [cm^2/sec]
    "Lead Pb+2"             : 0.945 * pow(10, -5),
    "Potassium K+"          : 1.960 * pow(10, -5),
    "Calcium Ca+2"          : 0.793 * pow(10, -5),
    "Sodium Na+"            : 1.330 * pow(10, -5),
    "Electrons in Silicon"  : 0.00025
}

d_L = 0.8 * pow(10, -4)
d_x = np.linspace(0, d_L, num=1000)
d_pos = 0.25 * np.sin(2 * np.pi * d_x / d_L) + 0.05 * np.sin(4 * np.pi * d_x / d_L)
d_neg = np.multiply(d_pos, -0.5)
d_potential_mat = np.vstack((d_pos, d_neg))
d_T = 1 / 600000
d_t_vec = np.array([0.6 * d_T, d_T])
debug_dict = {
    "ion_selection"         : ION_LIST[4],  # 0=Pb+2, 1=K+, 2=Ca+2, 3=Na+, 4=e-
    "ratchet_number"        : 2,
    "L"                     : d_L,
    "a1"                    : 0.25,
    "a2"                    : 0.05,
    "potential_profile"     : [d_L, d_x, d_t_vec, d_potential_mat],
    "flash_number"          : 1,
    "flash_mode"            : -0.5,
    "output_selection"      : 1
}


'''PHYSICAL CONSTANTS'''
ELECTRON_CHARGE = 1.6 * pow(10, -19)
NE = 1 * pow(10, 15)
BOLTZMANN_CONSTANT = 8.617333262 * pow(10, -5)
BOLTZMANN_CONSTANT_J = 1.380649 * pow(10, -23)
AVOGADRO_NUMBER = 6.0221409 * np.power(np.e, 23)
TEMPERATURE = 293


'''SIMULATION PARAMETERS'''
DATA_CSV_FILE = 'test_profiles 2.csv'
ALPHA = 0.5                             # amplitude factor for negative period of flashing ratchet
FLASHING_MODES = [0, -ALPHA]
BLANK_INT = 'blank'
NUMBER_OF_SIMULATIONS = 10000
SIGMA = (500 * pow(10, -4)) * (50 * pow(10, -7))
INTERVALS_FLASH_RATIO = 10
GRADIENT_COEFFICIENT = 0.0005
RESOLUTION = 1000
RATCHETS_IN_SYSTEM = 5
POINTS = 300
YELLOW = '#c49000'
PURPLE = '#892fba'

