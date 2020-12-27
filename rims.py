'''
rims.py

Simulator of the movement of ions in a ratchet based system

The program will simulate the movement of many ions with different initial conditions (IE) through the ratchet.
The IEs will be randomized, thus making the simulator Monte Carlo based.
Simulation will generate graphs plotting the location of the ions post simulation.  
The simulator works on every individual ion separately.
The ions location will be processed at time frames at a given rate.
Every such time frame is an iteration of the main loop, in which the location on the ion will be calculated and graphed.

Eran Weil
Itamar Eyal

Dr. Gideon Segev
Tel-Aviv university
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.constants import k
import progressbar
from datetime import datetime
import os
import random
from ion_parser import *
from ion_simulation import *
'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''
ELECREON_MODE = False
TEMPERATURE = 298           # T in kelvin
Kb = scipy.constants.k      # Boltzman constant
INTERVAL = 0.1

'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''

# input from user: ions in system, potential profile applied
# ions: atom isotop mols. example: Ag +1 2, H +1 1

print("-------------------------------------------------------")
print("     RIMS - Ratchet based Ion Movement Simulator")
print("-------------------------------------------------------")

input("\npress ENTER to begin...\n")

if ELECREON_MODE:
    print("In ELECTRON MODE, the simulated particles are elecreons.")

    ions_dict_list = [{
        'symbol': 'e',
        'chrage': periodictable.constants.electron_volt,
        'quantity': 1
        }]
    print(ions_dict_list)


else:

    print("-------------------------------------------------------")
    print("             Step 1- Configure the system")
    print("-------------------------------------------------------")

    print("Please specify the ions in the system and their quantities.")
    print("List ions in the following format: symbol charge mols\nSeperate ions by comma and press ENTER: \n")
    print("\nFor example: \tH +1 2, Ag +2 5\n")

    ions_verified = False

    while ions_verified == False:
        ion_str = input("Enter ions: ")
        if ion_str == '':
            continue
        ions_dict_list = ion_inputs_to_attributes(ion_str)
        ions_verified = verify_ions(ions_dict_list)

    print_ions(ions_dict_list)

print("-------------------------------------------------------")
print("             Step 2- Configure the ratchet")
print("-------------------------------------------------------")

print("Please describe the potential profile applied on the system.\n")
print("                    _          ")
print("        / |\         |         ")
print("      /   | \        | A[v]    ")
print("    /     |  \       |         ")
print("  /       |   \     _|         ")
print("                               ")
print("  \______/\____/               ")
print("    a[um]   b[um]            \n")

a = float(input("a[um] = "))
b = float(input("b[um] = "))
A = float(input("A[v] = "))

i = ion(1)
i.create_arena(a,b,A)

preview_potential = input("\nType y to see a preview of 1 period of the ratchet, Otherwise press ENTER to continue...")
if preview_potential == 'y':
    print("\nThe program will continue when you close the preview graph...")
    plt.show()
print("\nEnter ratchet flashing frequency in KHz:")
potential_frequency = int(input("Ratchet frequency [KHz] = "))




class rims:
    def __init__(self, ions_dict_list, potential_profile):
        self.ions = ions_dict_list
        self.potential_profile = potential_profile

        




