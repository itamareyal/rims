"""
defines.py

Hard coded data for the whole program
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
    #"Electron in Silicone"  : np.random.normal(2.5 * pow(10,5), 2 * pow(10,-4))
    "Electron in Silicone"  : 0.00025
    #"Electron in Silicone"  : 0.0036
}


'''PHYSICAL CONSTANTS'''
ELECTRON_CHARGE = 1.6 * pow(10, -19)
NE = 1 * pow (10, 21)
BOLTZMANN_CONSTANT = 8.617333262 * pow(10,-5)
BOLTZMANN_CONSTANT_J = 1.380649 * pow(10,-23)
AVOGADRO_NUMBER = 6.0221409 * pow(np.e,23)
TEMPERATURE = 293


'''SIMULATION PARAMETERS'''
ALPHA = 0.5                             # amplitude factor for negative period of flashing ratchet
FLASHING_MODES = [0, -ALPHA]
BLANK_INT = 'blank'
NUMBER_OF_SIMULATIONS = 10000
SIGMA = (500 * pow(10, -6)) * (50 * pow(10, -9))
INTERVALS_FLASH_RATIO = 10
GRADIENT_COEFFICIENT = 0.0005
RESOLUTION = 1000
RATCHETS_IN_SYSTEM = 5
POINTS = 300
YELLOW = '#c49000'
PURPLE = '#892fba'

