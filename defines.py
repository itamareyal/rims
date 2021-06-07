"""
defines.py

Hard coded data for the whole program
All units in cm / sec / Hz / kelvin / ev
"""

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import numpy as np
import csv


'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''

def load_settings(settings_file):
    settings_dict = {}
    with open(settings_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            settings_dict[row[0]] = row[1]
    return settings_dict

def detect_bool_from_str(string):
    if string.lower() in ['true', 'y', 'yes', 't']:
        return True
    return False


settings = load_settings('settings.csv')

'''PARTICLE PARAMETERS'''
ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electrons in Silicon"]

diffusion_coefficient_dict = {                      # units [cm^2/sec]
    # Taken from https://www.aqion.de/site/diffusion-coefficients
    "Lead Pb+2"             : 0.945 * pow(10, -5),
    "Potassium K+"          : 1.960 * pow(10, -5),
    "Calcium Ca+2"          : 0.793 * pow(10, -5),
    "Sodium Na+"            : 1.330 * pow(10, -5),
    # Taken from https://pubs.acs.org/doi/suppl/10.1021/acs.nanolett.7b03118/suppl_file/nl7b03118_si_001.pdf
    "Electrons in Silicon"  : 3 * pow(10, -4)
}
d_alpha = -1
d_dc = 0.8
d_f = 600_000
d_T = 1 / d_f
d_L = 0.8 * pow(10, -4)
d_a1 = 0.25
d_a2 = 0.05
d_x = np.linspace(start=0, stop=d_L, num=1000)
d_pos = d_a1 * np.sin(2 * np.pi * d_x / d_L) + d_a2 * np.sin(4 * np.pi * d_x / d_L)
d_neg = np.multiply(d_pos, d_alpha)
d_potential_mat = np.vstack((d_pos, d_neg))
d_t_vec = np.array([d_dc * d_T, d_T])

'''0=Pb+2, 1=K+, 2=Ca+2, 3=Na+, 4=e-'''
debug_dict = {
    "ion_selection"         : (ION_LIST[4], diffusion_coefficient_dict[ION_LIST[4]]),
    "ratchet_number"        : 2,
    "L"                     : d_L,
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
TEMPERATURE = 293


'''SIMULATION PARAMETERS'''
ENABLE_VIDEO = detect_bool_from_str(settings['ENABLE_VIDEO'])
ALPHA = -1
BLANK_INT = 'blank'
NUMBER_OF_SIMULATIONS = int(settings['NUMBER_OF_SIMULATIONS'])
STEADY_STATE_PERCENT_MARGIN = float(settings['STEADY_STATE_PERCENT_MARGIN'])
IONS_PER_THREAD = 100
NUMBER_OF_THREADS = int(NUMBER_OF_SIMULATIONS/IONS_PER_THREAD)
MAX_CYCLES = int(settings['MAX_CYCLES'])
SIGMA = (500 * pow(10, -4)) * (50 * pow(10, -7))
INTERVALS_FLASH_RATIO = 10
INTERVALS_FLASH_RATIO_ELECTRONS = 50
RESOLUTION = int(settings['RESOLUTION'])
RATCHETS_IN_SYSTEM = int(settings['RATCHETS_IN_SYSTEM'])
POINTS = 1000
MIN_NUM_SPEEDS_FOR_AVG = 10
YELLOW = '#c49000'
PURPLE = '#892fba'
BLUE = '#113b80'

