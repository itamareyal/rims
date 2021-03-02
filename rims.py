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
from video_gen import *
import numpy as np
import matplotlib.pyplot as plt
import sys


'''----------------------------------------------------------------------
                                DEFINES
----------------------------------------------------------------------'''


ION_LIST = ["Lead Pb+2", "Potassium K+", "Calcium Ca+2", "Sodium Na+", "Electron in Silicone"]
FLASHING_MODES = [0,-0.5]
BLANK_INT = 'blank'
NUMBER_OF_SIMULATIONS = 10000


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
        self.number_of_simulations = NUMBER_OF_SIMULATIONS
        self.start_time = datetime.now()
        self.path_for_output = 'RIMS output plots\\'+str(self.start_time.strftime("%x")).replace('/','-')+'_'+ str(self.start_time.strftime("%X")).replace(':','')+'\\'


    def create_ion(self):
        i = ion(self.ion, self.potential_profile_list, self.flash_frequency, self.flash_mode, self.dc, self.electric_field, self.potential_profile, self.path_for_output)
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
            x = np.linspace(0, L, num=RESOLUTION)
            V = a1 * np.sin(2*np.pi *x / L) + a2 * np.sin(4*np.pi *x / L)
            #dx = L/RESOLUTION


        else: #saw
            a = self.potential_profile_list[0]
            b = self.potential_profile_list[1]
            A = self.potential_profile_list[2]

            x = np.linspace(0,a+b, num=RESOLUTION)
            f1=A * np.divide(x,a)
            f2=A * np.divide((x-(a+b)),(-b))
            step = np.heaviside(x-a,1)
            V=  f1 -step*f1 + step* f2
            #dx = (a+b)/RESOLUTION
        E= -np.gradient(V,RESOLUTION * 0.00005)

        plt.figure(0)
        plt.plot(x,V, label="V(x) potential profile", color='#f5bc42')
        plt.plot(x,E, label=r"E(x) electric field = -$\nabla $V", color='#892fba')
        plt.suptitle('RIMS: Ratchet potential profile', fontsize=14, fontweight='bold')
        plt.xlabel(r"X [$\mu $m]")
        plt.ylabel(r"V [v] , E [v/$\mu $m]")
        plt.legend(loc="upper left")
        self.save_plots('Ratchet potential profile',0)
        self.electric_field = E
        self.potential_profile = V
        return


    def save_plots(self, name,fig_number):
        if not os.path.exists(self.path_for_output):
            os.makedirs(self.path_for_output)
        plt.figure(fig_number)
        plt.savefig(self.path_for_output+name+'.jpeg')
        print(name+' saved to output plots.')


    def create_trace_file(self, ion_subject):
        with open(self.path_for_output + 'simulation trace.csv',newline='', mode='a') as csv_file:
            writer = csv.writer(csv_file,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['X0[um]', 'X'+str(ion_subject.points)+'[um]'])


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
        print("Plotting results...\nClose the plots to continue")
        plt.figure(1)
        weights = np.ones_like(x_results)/float(len(x_results))
        plt.hist(x_results, weights=weights, bins=RESOLUTION, label=r"X [$\mu $m]")
        # mn, mx = plt.xlim()
        # plt.xlim(mn, mx)
        # kde_xs = np.linspace(mn, mx, 300)
        # kde = st.gaussian_kde(x_results)
        # plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
        #plt.legend(loc="upper left")
        plt.ylabel('Density')
        plt.xlabel(r'X [$\mu $m]')
        plt.title(r"Histogram of distribution x axis: $\rho $(x)")
        self.save_plots('Distribution x axis histogram',1)
        plt.show()


    def create_video(self):
        print("\nRIMS simulation in progress...")
        print("\nCreating frames for video...")
        frames = 100
        ion_subject = self.create_ion()
        for frame in range(frames):
            frame+=1
            x_results = []
            simulation_count = 0

            if not os.path.exists(self.path_for_output + 'frames'):
                os.makedirs(self.path_for_output + 'frames')

            self.create_trace_file(ion_subject)

            while simulation_count < self.number_of_simulations:
                ion_subject = self.create_ion()
                ion_subject.points = frame
                x_results.append(ion_subject.simulate_ion())
                simulation_count += 1

                prog = simulation_count * 100 / int(self.number_of_simulations)
                sys.stdout.write("\r%d%% " % prog)
                sys.stdout.flush()

                self.write_to_trace_file(ion_subject)
            plt.figure(frame)
            weights = np.ones_like(x_results) / float(len(x_results))
            plt.hist(x_results, weights=weights, bins=RESOLUTION, label=r"X [$\mu $m]", range=(-5,5))
            # mn, mx = plt.xlim()
            # plt.xlim(mn, mx)
            # kde_xs = np.linspace(mn, mx, 301)
            # kde = st.gaussian_kde(x_results)
            # plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
            #plt.legend(loc="upper left")
            plt.ylim(0,0.035)
            plt.ylabel('Density')
            plt.xlabel(r'X [$\mu $m]')
            plt.title(r"Histogram of distribution x axis: $\rho $(x,t)")
            plt.suptitle('t = '+str(ion_subject.interval * ion_subject.intervals_count)[0:8]+r' [$\mu $sec]', fontsize=10)
            self.save_plots('frames\\frame_'+str(frame), frame)
            plt.close(frame)
        print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
        self.create_log_file(ion_subject)
        print("Simulation log file created and saved.\n")



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
        f.write("\tintervals (delta_t): "+str(ion_subject.interval)+"e-6[sec]\n")
        f.write("\tfriction coefficient (gamma): " + str(ion_subject.gamma) + "\n")
        f.write("\tresolution: " + str(RESOLUTION) + "\n")

        f.close()

def input_check_int(msg, range):
    val = BLANK_INT
    while val not in range:
        try:
            val = int(input(msg))
        except ValueError:
            print("\tPlease enter an integer as specified above")
            continue
    return val

def input_check_float(msg):
    val = BLANK_INT
    clear = False
    while not clear:
        try:
            val = float(input(msg))
            clear = True
        except ValueError:
            print("\tPlease enter an integer or a float as specified above")
            continue
    return val

'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''


print("-------------------------------------------------------")
print("     RIMS - Ratchet based Ion Movement Simulator")
print("-------------------------------------------------------")

input("\npress ENTER to begin...\n")
def execution():
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
    number_selection = input_check_int("Enter your selection:", range(1,7))
    # while number_selection not in range(7):
    #     try:
    #         number_selection = int(input("Enter your selection:")) -1
    #     except ValueError:
    #         continue

    if number_selection == 6:
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
        output_selection =1

    else:

        ion_selection = ION_LIST[number_selection-1]

        print("\nIon selected: "+ion_selection)

        print("-------------------------------------------------------")
        print("             Step 2- Configure the ratchet")
        print("-------------------------------------------------------")

        print("\nEnter ratchet function:\n\t1)Saw wave\n\t2)Double Sin\n")
        ratchet_number = input_check_int("Ratchet function number = ", [1,2])
        # ratchet_number = BLANK_INT
        # while ratchet_number not in [1,2]:
        #     try:
        #         ratchet_number = int(input("Ratchet function = "))
        #     except ValueError:
        #         continue
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

            #a = float(input("\ta[um] = "))
            a = input_check_float("\ta[um] = ")
            b = input_check_float("\tb[um] = ")
            A = input_check_float("\tA[v] = ")
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
            L = input_check_float("\tL[um] = ")
            a1 = input_check_float("\ta1[v] = ")
            a2 = input_check_float("\ta2[v] = ")
            potential_profile = [L, a1, a2, ratchet_number]


        print("\nEnter ratchet flashing frequency in Hz: (can use factors of K,M,G)")
        flash_frequency = -1
        tries =0
        while flash_frequency <= 0:
            if tries>0:
                print("\tfrequency takes values equal or larger than 1")
            try:
                raw_input = input("Ratchet frequency [Hz] = ")
                converted_to_int = raw_input
                converted_to_int = converted_to_int.replace('k','').replace('K','').replace('M','').replace('G','')
                flash_frequency = int(converted_to_int)
                if 'k' in raw_input or 'K' in raw_input:
                    flash_frequency *= 1000
                if 'M' in raw_input:
                    flash_frequency *= 1000000
                if 'G' in raw_input:
                    flash_frequency *= 1000000000

            except ValueError:
                print("\tPlease enter an integer as specified above")
                tries += 1
                continue
            tries+=1

        print("\nEnter ratchet duty cycle from 0-1:")
        dc = -1
        tries =0
        while dc < 0 or dc > 1:
            if tries>0:
                print("\tdc takes float values from (0-1)")
            try:
                dc = float(input("DC = "))
            except ValueError:
                tries += 1
                continue
            tries+=1

        print("\nEnter flashing mode number:\n\t1)ON/OFF\n\t2)+/-\n")
        flash_number = input_check_int("Flashing mode = ",[1,2])
        flash_mode = FLASHING_MODES[flash_number-1]

        print("-------------------------------------------------------")
        print("             Step 3- Outputs selection")
        print("-------------------------------------------------------")
        print("\nEnter desired output combination:\n\t1)Histogram (about 30sec to generate)\n\t2)Video (about 40min to generate)")
        output_selection = input_check_int("Enter your selection:",[1,2])

        print("\n-------------------------------------------------------")
        print("             Starting simulation")
        print("-------------------------------------------------------")
        print("\nBuilding Monte Carlo calculation system")
        try:
            print("Initial location of every ion is uniformly randomized over [0,"+str(L)+"um]\n")
        except NameError:
            print("Initial location of every ion is uniformly randomized over [0," + str(a+b) + "um]\n")

    r= rims(ion_selection, potential_profile, flash_frequency, flash_mode, dc)
    r.electric_field()
    if output_selection==1:
        r.create_histogram()
    elif output_selection==2:
        r.create_video()
        generate_video_from_frames(r.path_for_output + 'frames', 'density over time.avi')

    print("\n-------------------------------------------------------")
    print("                 Simulation over")
    print("-------------------------------------------------------")
    rerun = input("\nPress y and ENTER to run an new simulation, otherwise press ENTER...\n")
    if rerun == 'y':
        return 1
    return 0

while execution():
    continue


        




