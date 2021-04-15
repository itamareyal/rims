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
ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electron in Silicone"]

diffusion_coefficient_dict = {              # Deff dictionary. [m^2/sec]
    "Lead Pb+2"             : 0.945,
    "Potassium K+"          : 1.960,
    "Calcium Ca+2"          : 0.793,
    "Sodium Na+"            : 1.330,
    "Electron in Silicone"  : 0.00025       # cm^2 / sec
}

debug_dict = {
    "ion_selection"         : ION_LIST[4],  # 0=Pb+2, 1=K+, 2=Ca+2, 3=Na+, 4=e-
    "ratchet_number"        : 2,
    "L"                     : 0.8 * pow(10, -4),
    "a1"                    : 0.25,
    "a2"                    : 0.05,
    "potential_profile"     : [0.8 * pow(10, -4), 0.25, 0.05, 2],
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

