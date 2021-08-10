import numpy as np
import csv

"""
defines.py

Hard coded data for the whole program
All units in cm / sec / Hz / kelvin / ev
"""

'''----------------------------------------------------------------------
                        EXTRACTION OF SETTINGS
----------------------------------------------------------------------'''


def load_settings(settings_file):
    """
    Loads data from csv file into dictionary
    """
    settings_dict = {}
    with open(settings_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0][0] in ['#', '\n']:
                continue
            value = row[1]
            try:
                value = float(row[1])
            except ValueError:
                pass
            settings_dict[row[0]] = value
    return settings_dict

def detect_bool_from_str(string):
    if string.lower() in ['true', 'y', 'yes', 't']:
        return True
    return False


settings = load_settings('settings.csv')
diffusion_coefficient_dict = load_settings('diffusion_coefficients.csv')

'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''

'''VERSION'''
VERSION = '1.9.5'

'''DEBUG PARAMETERS'''
d_alpha = -1
d_dc = 0.6
d_f = 900_000
d_T = 1 / d_f
d_L = 0.8 * pow(10, -4)
d_a1 = 0.25
d_a2 = 0.05
d_x = np.linspace(start=0, stop=d_L, num=1000)
d_pos = d_a1 * np.sin(2 * np.pi * d_x / d_L) + d_a2 * np.sin(4 * np.pi * d_x / d_L)
d_neg = np.multiply(d_pos, d_alpha)
d_potential_mat = np.vstack((d_pos, d_neg))
d_t_vec = np.array([d_dc * d_T, d_T])
d_potential_profile = [d_L, d_x, d_t_vec, d_potential_mat]
d_ion = 'Electrons in Silicon'
d_diffusion = diffusion_coefficient_dict[d_ion]
d_ion_selection = (d_ion, d_diffusion)

'''PHYSICAL CONSTANTS'''
ELECTRON_CHARGE = 1.6 * pow(10, -19)
BOLTZMANN_CONSTANT = 8.617333262 * pow(10, -5)
BOLTZMANN_CONSTANT_J = 1.380649 * pow(10, -23)
TEMPERATURE = 293

'''SIMULATION PARAMETERS'''
ENABLE_VIDEO = detect_bool_from_str(settings['ENABLE_VIDEO'])
CREATE_TRACE = detect_bool_from_str(settings['CREATE_TRACE'])
ALPHA = -1
PARTICLES_SIMULATED = int(settings['PARTICLES_SIMULATED'])
STEADY_STATE_PERCENT_MARGIN = float(settings['STEADY_STATE_PERCENT_MARGIN'])
IONS_PER_THREAD = 100
NUMBER_OF_THREADS = int(PARTICLES_SIMULATED/IONS_PER_THREAD)
MAX_CYCLES = int(settings['MAX_CYCLES'])
INTERVALS_FLASH_RATIO = 10
INTERVALS_FLASH_RATIO_ELECTRONS = 50
RESOLUTION = int(settings['RESOLUTION'])
RATCHETS_IN_SYSTEM = int(settings['RATCHETS_IN_SYSTEM'])
MIN_MEASUREMENTS_FOR_SS = 10
OVERWRITE_DELTA_T = detect_bool_from_str(settings['OVERWRITE_DELTA_T'])
DELTA_T = float(settings['DELTA_T']) * pow(10, -6)

'''KEDEM PROBLEM PARAMETERS'''
NE = 1 * pow(10, 15)
SIGMA = (500 * pow(10, -4)) * (50 * pow(10, -7))

'''COLORS'''
YELLOW = '#c49000'
PURPLE = '#892fba'
BLUE = '#113b80'

