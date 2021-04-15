""""
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
"""

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import sys
from ion_simulation import *
from current_calc import *
from video_gen import *
from output import *
from defines import *


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''


class rims:
    def __init__(self, ion_subject, potential_profile, flash_frequency, flash_mode, dc):
        """
        Instance holding attributes of the whole simulation
        :param ion_subject: ion type to be simulated
        :param potential_profile: list holding potential profile shape [L,a1,a2]
        :param flash_frequency: flashing rate of the ratchet
        :param flash_mode: determines the profile when in 1-dc. Either completely off or negative and amp by ALPHA
        :param dc: duty cycle of positive part of the potential profile
        """
        self.ion = ion_subject
        self.current = 0
        self.electric_field = 0
        self.potential_profile = 0
        self.potential_profile_list = potential_profile
        self.flash_frequency = flash_frequency
        self.flash_mode = flash_mode
        self.dc = dc
        self.number_of_simulations = NUMBER_OF_SIMULATIONS
        self.start_time = datetime.now()
        self.path_for_output = r'RIMS output plots/' + get_time_stamp(self) + r'/'

    def create_ion(self):
        """
        Creates ion class instance to be simulated
        :return: ion instance
        """
        i = ion(self.ion, self.potential_profile_list, self.flash_frequency, self.flash_mode, self.dc, self.electric_field, self.potential_profile, self.path_for_output) #Add self.temperature if necessary
        i.create_arena()
        i.get_intervals()
        i.get_gamma()
        return i

    def get_electric_field(self):
        """
        Derives the electric field from the potential, E(x).
        saves E(x) as attribute self.electric field and V(x) as self.potential_profile.
        """
        if self.potential_profile_list[3] == 2:     # sin
            L = self.potential_profile_list[0]      # profile length x axis [um]
            a1 = self.potential_profile_list[1]     # amp of first sin wave [V]
            a2 = self.potential_profile_list[2]     # amp of second sin wave [V]
            x = np.linspace(0, L, num=RESOLUTION)
            V = a1 * np.sin(2*np.pi *x / L) + a2 * np.sin(4*np.pi *x / L)

        else:                                       # saw
            a = self.potential_profile_list[0]      # narrow left part of saw [um]
            b = self.potential_profile_list[1]      # thick right part of saw [um]
            A = self.potential_profile_list[2]      # amp of saw wave [V]
            x = np.linspace(0, a+b, num=RESOLUTION)
            f1 = A * np.divide(x, a)
            f2 = A * np.divide((x-(a+b)), (-b))
            step = np.heaviside(x-a, 1)
            V = f1 - step * f1 + step * f2

        E = -np.gradient(V, 0.001)
        self.electric_field = E
        self.potential_profile = V

        '''plot E & V'''
        plot_potential_profile(self, x, V, E)
        return

    def create_histogram(self):
        """
        Running ions in the system and collecting thier location after the simulation.
        Creating a histogram from all collected locations
        """
        print("\nRIMS simulation in progress...")

        x_results = []
        simulation_count = 0
        ion_subject = self.create_ion()
        create_trace_file(self, ion_subject)
        vt_list = []
        # hops_list = []
        '''main simulation loop'''
        while simulation_count < self.number_of_simulations:
            ion_subject = self.create_ion()
            x_results.append(ion_subject.simulate_ion())
            vt_list.append(ion_subject.vi)

            simulation_count += 1
            progress = simulation_count * 100 / int(self.number_of_simulations)
            sys.stdout.write("\r%d%%" % progress)
            sys.stdout.flush()
            write_to_trace_file(self, ion_subject)

        '''calculation of particles velocity and current at steady state'''
        vt_vector = np.array(vt_list)
        vt_av = np.average(vt_vector)
        vt_over_T = vt_av * ion_subject.interval * self.flash_frequency
        self.current = get_current(-vt_over_T, NE, SIGMA, ELECTRON_CHARGE)

        print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
        create_log_file(self, ion_subject)
        print("Simulation log file created and saved.\n")

        '''plotting distribution of particles histogram'''
        plot_id = create_unique_id()
        plt.figure(plot_id)
        weights = np.ones_like(x_results) / float(len(x_results))
        plt.hist(x_results, weights=weights, bins=RESOLUTION, label=r"X [$\mu $m]")
        plt.ylabel('Density')
        plt.xlabel(r'X [$\mu $m]')
        plt.title(r"Histogram of distribution x axis: $\rho $(x)")
        save_plots(self, 'Distribution x axis histogram', plot_id)
        return

    def create_video(self):
        print("\nRIMS simulation in progress...")
        print("\nCreating frames for video...")
        frames = 100
        ion_subject = self.create_ion()
        for frame in range(frames):
            frame += 1
            x_results = []
            simulation_count = 0

            if not os.path.exists(self.path_for_output + 'frames'):
                os.makedirs(self.path_for_output + 'frames')

            create_trace_file(self, ion_subject)

            while simulation_count < self.number_of_simulations:
                ion_subject = self.create_ion()
                ion_subject.points = frame
                x_results.append(ion_subject.simulate_ion())

                simulation_count += 1

                prog = simulation_count * 100 / int(self.number_of_simulations)
                sys.stdout.write("\r%d%% " % prog)
                sys.stdout.flush()

                write_to_trace_file(self, ion_subject)

            plt.figure(frame)
            weights = np.ones_like(x_results) / float(len(x_results))
            plt.hist(x_results, weights=weights, bins=RESOLUTION, label=r"X [$\mu $m]", range=(-5, 5))
            plt.ylim(0, 0.035)
            plt.ylabel('Density')
            plt.xlabel(r'X [$\mu $m]')
            plt.title(r"Histogram of distribution x axis: $\rho $(x,t)")
            plt.suptitle('t = ' + str(ion_subject.interval * ion_subject.intervals_count)[0:8] + r' [$\mu $sec]',
                         fontsize=10)
            save_plots(self, 'frames\\frame_' + str(frame), frame)
            plt.close(frame)
        print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
        create_log_file(self, ion_subject)
        print("Simulation log file created and saved.\n")


def input_check_int(msg, desired_range):
    val = BLANK_INT
    while val not in desired_range:
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


def get_time_stamp(rims_object):
    time_start = rims_object.start_time
    ts = str(time_start.strftime("%x")).replace('/', '-')+'_' \
        + str(rims_object.start_time.strftime("%X")).replace(':', '')
    return ts


'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''

print("-------------------------------------------------------")
print("     RIMS - Ratchet based Ion Movement Simulator")
print("-------------------------------------------------------")

input("\npress ENTER to begin...\n")


def execution(dc_sample, f_sample):
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
    # number_selection = input_check_int("Enter your selection:", range(1,7))
    number_selection = 6

    if number_selection == 6:
        ion_selection = ION_LIST[4]
        ratchet_number = 2
        L = 0.8 * pow(10, -4)
        a1 = 0.25
        a2 = 0.05
        potential_profile = [L, a1, a2, ratchet_number]
        flash_frequency = f_sample
        dc = dc_sample
        flash_number = 1
        flash_mode = -0.5
        output_selection = 1

    else:

        ion_selection = ION_LIST[number_selection-1]

        print("\nIon selected: "+ion_selection)

        print("-------------------------------------------------------")
        print("             Step 2- Configure the ratchet")
        print("-------------------------------------------------------")

        print("\nEnter ratchet function:\n\t1)Saw wave\n\t2)Double Sin\n")
        ratchet_number = input_check_int("Ratchet function number = ", [1,2])
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

            a = input_check_float("\ta[um] = ")
            b = input_check_float("\tb[um] = ")
            L = a + b
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
        tries = 0
        while flash_frequency <= 0:
            if tries>0:
                print("\tfrequency takes values equal or larger than 1")
            try:
                raw_input = input("Ratchet frequency [Hz] = ")
                converted_to_int = raw_input
                converted_to_int = converted_to_int.replace('k', '').replace('K', '').replace('M', '').replace('G', '')
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
            tries += 1

        print("\nEnter ratchet duty cycle from 0-1:")
        dc = -1
        tries = 0
        while dc < 0 or dc > 1:
            if tries > 0:
                print("\tdc takes float values from (0-1)")
            try:
                dc = float(input("DC = "))
            except ValueError:
                tries += 1
                continue
            tries += 1

        print("\nEnter flashing mode number:\n\t1)ON/OFF\n\t2)+/-\n")
        flash_number = input_check_int("Flashing mode = ", [1, 2])
        flash_mode = FLASHING_MODES[flash_number-1]

        print("-------------------------------------------------------")
        print("             Step 3- Outputs selection")
        print("-------------------------------------------------------")
        print("\nEnter desired output combination:\n\t1)Histogram (about 30sec to generate)\n\t"
              "2)Video (about 40min to generate)")
        output_selection = input_check_int("Enter your selection:", [1, 2])

        print("\n-------------------------------------------------------")
        print("             Starting simulation")
        print("-------------------------------------------------------")
        print("\nBuilding Monte Carlo calculation system")

        print("Initial location of every ion is uniformly randomized over [0," + str(L) + "cm]\n")

    r = rims(ion_selection, potential_profile, flash_frequency, flash_mode, dc)
    r.get_electric_field()
    if output_selection == 1:
        r.create_histogram()
    elif output_selection == 2:
        r.create_video()
        generate_video_from_frames(r.path_for_output + 'frames', 'density over time.avi')

    print_log_file(r)

    print("\n-------------------------------------------------------")
    print("                 Simulation over")
    print("-------------------------------------------------------")
    # rerun = input("\nPress y and ENTER to run an new simulation, otherwise press ENTER...\n")
    # # if rerun == 'y':
    # #     return 1
    return r.current


def create_i_of_dc():
    """
    Runs 20 different dc samples to create I(dc) graph
    """
    dc_list = []
    current_list = []
    plot_uid = create_unique_id()
    runs = 40
    for dc_sample in range(0, runs):
        current_calculated = execution(dc_sample/runs, 600000)
        dc_list.append(dc_sample/runs)
        current_list.append(current_calculated)
    plt.figure(plot_uid)
    plt.plot(dc_list, current_list, label="I(duty cycle)", color='#f5bc42')
    plt.suptitle('RIMS: current changing over DC at 600KHz', fontsize=12, fontweight='bold')
    plt.xlabel(r"DC")
    plt.ylabel(r"I [AMP]")

    max_current = max(current_list)
    max_i = current_list.index(max_current)
    plt.plot(dc_list[max_i], max_current, 'g^')
    plt.text(dc_list[max_i], max_current, str(max_current))

    min_current = min(current_list)
    min_i = current_list.index(min_current)
    plt.plot(dc_list[min_i], min_current, 'rv')
    plt.text(dc_list[min_i], min_current, str(min_current))

    plt.axhline(color='r')
    plt.savefig('i_dc ' + plot_uid + '.jpeg')
    print('I_dc saved to output plots.')


def create_heat_map():
    resolution_f = 10
    resolution_dc = 10
    dc_list = [dc/resolution_dc for dc in range(resolution_dc, 0, -1)]
    f_list = [f * 1000 for f in range(100, 1100, 100)]

    plot_uid = create_unique_id()
    matrix = np.zeros(shape=(resolution_f, resolution_dc))
    for i_f in range(0, resolution_f):
        current_vector = np.zeros(resolution_dc)
        for i_dc in range(0, resolution_dc):
            dc_input = dc_list[i_dc]
            f_input = f_list[i_f]
            current_calculated = execution(dc_input, f_input)
            current_vector[i_dc] = current_calculated
        matrix[resolution_f - i_f - 1] = current_vector

    fig, ax = plt.subplots()

    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    im, _ = heatmap(matrix.transpose(), dc_list, f_list_label, ax=ax, vmin=0,
                    cmap="PuOr", cbarlabel="I(DC, frequency)")
    plt.tight_layout()
    ax.set_title('RIMS: I(DC, frequency) heat-map')
    ax.set_ylabel('Duty Cycle')
    ax.set_xlabel('Ratchet Frequency')

    file_name = 'i_dc_f_heat-map ' + plot_uid + '.jpeg'
    plt.savefig(file_name)
    print('I_dc_f_heat map saved to main directory as ' + file_name)
    plt.close(plot_uid)
    return


create_i_of_dc()
# create_heat_map()
