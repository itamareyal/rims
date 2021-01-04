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


import scipy.stats as st
from datetime import datetime
from ion_simulation import *
import numpy as np
import matplotlib.pyplot as plt
import sys


'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''


ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electron in Silicone"]
FLASHING_MODES = [0,-0.5]


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''


class rims:
    def __init__(self, ion, potential_profile, flash_frequency, flash_mode, dc):
        self.ion = ion
        self.potential_profile_list = potential_profile
        self.flash_frequency = flash_frequency
        self.flash_mode = flash_mode
        self.dc = dc
        self.number_of_simulations = 10000
        self.start_time = datetime.now()
        self.path_for_output = 'RIMS output plots\\'+str(self.start_time.strftime("%x")).replace('/','-')+'_'+ str(self.start_time.strftime("%X")).replace(':','')+'\\'


    def create_ion(self):
        i = ion(self.ion, self.potential_profile_list, self.flash_frequency, flash_mode, self.dc, self.electric_field, self.potential_profile, self.path_for_output)
        i.create_arena()
        i.get_intervals()
        i.get_gamma()
        return i


    def electric_field(self):
        # Description: derives the electric field from the potential, E(x).
        # Parameters: self
        # Return: saves E(x) as attribute self.electric field and V(x) as self.potential_profile.

        if self.potential_profile_list[3] == 2: #sin
            L = self.potential_profile_list[0]
            a1 = self.potential_profile_list[1]
            a2 = self.potential_profile_list[2]
            x = np.linspace(0, L)
            V = a1 * np.sin(2*np.pi *x / L) + a2 * np.sin(4*np.pi *x / L)
            E = -np.gradient(V)

        else: #saw
            a = self.potential_profile_list[0]
            b = self.potential_profile_list[1]
            A = self.potential_profile_list[2]

            x = np.linspace(0,a+b)
            f1=A * np.divide(x,a)
            f2=A * np.divide((x-(a+b)),(-b))
            step = np.heaviside(x-a,1)
            V=  f1 -step*f1 + step* f2
            E= -np.gradient(V)

        plt.figure(0)
        plt.plot(x,V, label="V(x)")
        plt.plot(x,E, label="E(x)")
        plt.suptitle('RIMS: Ratchet potential profile', fontsize=14, fontweight='bold')
        plt.xlabel("X [um]")
        plt.ylabel("V [v]")
        plt.legend(loc="upper left")
        self.save_plots('Ratchet potential profile',0)

        self.electric_field = E
        self.potential_profile = V
        return


    def save_plots(self, name,fig_number):
        if not os.path.exists(self.path_for_output):
            os.makedirs(self.path_for_output)
        plt.figure(fig_number)
        plt.savefig(self.path_for_output+name+'.png')
        print(name+' saved to output plots.')


    def create_trace_file(self, ion_subject):
        with open(self.path_for_output + 'simulation trace.csv',newline='', mode='a') as csv_file:
            writer = csv.writer(csv_file,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['X0', 'X'+str(ion_subject.points)])


    def write_to_trace_file(self, ion_subject):
        with open(self.path_for_output + 'simulation trace.csv', newline='', mode='a') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            # writer.writerow(['X0', 'X'+str(self.points)])
            writer.writerow([ion_subject.x0, ion_subject.L * ion_subject.arena_count + ion_subject.loc])


    def create_histogram(self):
        print("\nRIMS simulation in progress...")

        x_results = []
        simulation_count =0
        ion_subject = self.create_ion()

        self.create_trace_file(ion_subject)

        while simulation_count < self.number_of_simulations:
            ion_subject = self.create_ion()
            x_results.append(ion_subject.simulate_ion())
            simulation_count += 1
            prog = simulation_count*100/int(self.number_of_simulations)
            sys.stdout.write("\r%d%%" % prog)
            sys.stdout.flush()
            self.write_to_trace_file(ion_subject)


        print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
        self.create_log_file(ion_subject)
        print("Simulation log file created and saved.\n")
        print("Plotting results...\n")
        plt.figure(1)
        plt.hist(x_results, density=True, bins=30, label="X[um]")
        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        kde_xs = np.linspace(mn, mx, 301)
        kde = st.gaussian_kde(x_results)
        plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
        plt.legend(loc="upper left")
        plt.ylabel('Probability')
        plt.xlabel('X[um]')
        plt.title("Histogram")
        self.save_plots('Distribution x axis histogram',1)
        plt.show()


    def create_log_file(self, ion_subject):
        if not os.path.exists(self.path_for_output):
            os.makedirs(self.path_for_output)
        f = open(self.path_for_output+"RIMS simulation log.txt", "a")
        f.write("RIMS simulation log\n\n")
        f.write("\ttime created: "+str(self.start_time)+"\n")
        f.write("\ttest duration: "+str(datetime.now() - self.start_time)+"\n")

        f.write("\n\tparticles in the system: " + self.ion + "\n")
        f.write("\tdiffusion coefficient: " + str(ion_subject.diffusion) + "[m^2/cm] /10^-9\n")

        f.write("\nRatchet potential profile\n")
        if self.potential_profile_list[3]==2: #sin
            f.write("\tfunction: double sin wave \n")
            f.write("\tamplitude: " + str( self.potential_profile_list[1] +self.potential_profile_list[2]) + "[V]\n")

        else: #saw
            f.write("\tfunction: saw wave \n")
            f.write("\tamplitude: " + str(self.potential_profile_list[1]) + "[V]\n")
        f.write("\twidth: " + str(ion_subject.L) + "[um]\n")
        f.write("\tfrequency: " + str(ion_subject.flash_frequency) + "[Hz]\n")
        f.write("\tperiod: " + str(ion_subject.flash_period) + "[sec]\n")
        f.write("\tduty cycle: " + str(self.dc) +"\n")
        if self.flash_mode == 0: #ON/OFF
            f.write("\tflash mode: ON/OFF\n")
        else:
            f.write("\tflash mode: + / -\n")

        f.write("\nSimulation settings\n")
        f.write("\tparticles simulated: " + str(self.number_of_simulations) + "\n")
        f.write("\tmeasurements per particle: " + str(ion_subject.points) + "\n")
        f.write("\tintervals (delta_t): "+str(ion_subject.interval)+"[sec]\n")
        f.write("\tfriction coefficient (gamma): " + str(ion_subject.gamma) + "\n")

        f.close()


'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''


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
print("\t6) debug")
number_selection = int(input("Enter your selection:")) -1

if number_selection == 5:
    ion_selection = ION_LIST[4]
    ratchet_number =2
    L = 0.8
    a1 = 0.25
    a2 = 0.05
    potential_profile = [L, a1, a2, ratchet_number]
    flash_frequency = 600000
    dc = 0.6
    flash_number = 1
    flash_mode = -0.5

else:

    ion_selection = ION_LIST[number_selection]

    print("\nIon selected: "+ion_selection)

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


    print("\nEnter ratchet flashing frequency in Hz:")
    flash_frequency = int(input("Ratchet frequency [Hz] = "))
    print("\nEnter ratchet duty cycle from 0-1:")
    dc = float(input("DC = "))
    print("\nEnter flashing mode number:\n\t1)ON/OFF\n\t2)+/-\n")
    flash_number = int(input("Flashing mode = ")) -1
    flash_mode = FLASHING_MODES[flash_number]


r= rims(ion_selection, potential_profile, flash_frequency, flash_mode, dc)
r.electric_field()
r.create_histogram()

        




