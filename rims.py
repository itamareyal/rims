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
import scipy.stats as st
import progressbar
from datetime import datetime
import os
import random
#from ion_parser import *
from ion_simulation import *

'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''

ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electron in Silicone"]
FLASHING_MODES = [0,-1]


'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''

# input from user: ions in system, potential profile applied
# ions: atom isotop mols. example: Ag +1 2, H +1 1

print("-------------------------------------------------------")
print("     RIMS - Ratchet based Ion Movement Simulator")
print("-------------------------------------------------------")

input("\npress ENTER to begin...\n")

print("-------------------------------------------------------")
print("             Step 1- Configure the system")
print("-------------------------------------------------------")

print("\nSelect an ion to be simulated from the following list:")
print("\t1) Lead Pb+2")
print("\t2) Potassium K+")
print("\t3) Calcium Ca+2")
print("\t4) Sodium Na+")
print("\t5) Electron in Silicone")
number_selection = int(input("Enter your selection:")) -1
ion_selection = ION_LIST[number_selection]

print("\nIon selected: "+ion_selection)
input("Press Enter to confirm...")

print("-------------------------------------------------------")
print("             Step 2- Configure the ratchet")
print("-------------------------------------------------------")

print("\nEnter ratchet function:\n\t1)Saw wave\n\t2)Double Sin\n")
ratchet_number = int(input("Ratchet function = "))
if ratchet_number ==1:

    print("Please describe the potential profile applied on the system.\n")
    print("                    _          ")
    print("        / |\         |         ")
    print("      /   | \        | A[v]    ")
    print("    /     |  \       |         ")
    print("  /       |   \     _|         ")
    print("                               ")
    print("  \______/\____/               ")
    print("    a[um]   b[um]            \n")

    a = float(input("\ta[um] = "))
    b = float(input("\tb[um] = "))
    A = float(input("\tA[v] = "))
    potential_profile = [a, b, A, ratchet_number]

else:
    print("Please enter ratchet sin wave parameters.\n")
    print("                    _          ")
    print("        / |\        _| a2[v]   ")
    print("      /   | \        |         ")
    print("    /     |  \       | a1[v]   ")
    print("  /       |   \     _|         ")
    print("                               ")
    print("  \____________/               ")
    print("        L[um]                \n")
    print("qV(x) = a1 * sin(2pi * x / L) + a2 * sin(4pi * x / L)\n")
    L = float(input("\tL[um] = "))
    a1 = float(input("\ta1[v] = "))
    a2 = float(input("\ta2[v] = "))
    potential_profile = [L, a1, a2, ratchet_number]


print("\nEnter ratchet flashing frequency in KHz:")
flash_frequency = int(input("Ratchet frequency [KHz] = ")) *1000
print("\nEnter flashing mode number:\n\t1)ON/OFF\n\t2)+/-\n")
flash_number = int(input("Flashing mode = ")) -1
flash_mode = FLASHING_MODES[flash_number]


class rims:
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode):
        self.ion = ion
        self.potential_profile = potential_profile
        self.flash_frequency = flash_frequency
        self.flash_mode = flash_mode
        self.number_of_simulations = 100

    def create_ion(self):
        self.start_time = datetime.now()
        print(self.ion)
        i = ion(self.ion, self.potential_profile, self.flash_frequency, flash_mode)
        i.create_arena()
        i.get_intervals()
        i.get_gamma()
        print("\nIon profile for simulation created.")
        return i

    def create_histogram(self, ion):
        print("\nRIMS simulation in progress...")

        bar = progressbar.ProgressBar(maxval=self.number_of_simulations,
        widgets=[progressbar.Bar('=', '\t[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        x_results = []
        simulation_count =0
        while simulation_count < self.number_of_simulations:
            x_results.append(ion.simulate_ion())
            simulation_count += 1
            bar.update(simulation_count)
        print("\nDone!")
        print("\nPlotting results...")
        plt.hist(x_results, density=True, bins=30, label="Data")
        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        kde_xs = np.linspace(mn, mx, 301)
        kde = st.gaussian_kde(x_results)
        plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
        plt.legend(loc="upper left")
        plt.ylabel('Probability')
        plt.xlabel('Data')
        plt.title("Histogram")
        print("\ndone.")
        plt.show()


r= rims(ion_selection, potential_profile, flash_frequency, flash_mode)

r.create_histogram(r.create_ion())

        




