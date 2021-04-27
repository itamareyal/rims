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


from ion_simulation import *
from current_calc import *
from video_gen import *
from output import *
from defines import *
from interface import *


'''----------------------------------------------------------------------
                         CLASS IMPLEMENTATION
----------------------------------------------------------------------'''


class rims:
    def __init__(self, ion_subject, potential_profile, flash_frequency, flash_mode, dc, fast_mode):
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
        self.fast_mode = fast_mode
        self.steady_state_matrix = []
        self.steady_state = False

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
        if not self.fast_mode:
            plot_potential_profile(self, x, V, E)
        return

    def check_for_steady_state(self, vt_list):
        if len(vt_list) < 5:
            return
        last_vi_margin = vt_list[-1] / 10
        for i in range(2, 6, 1):
            if (vt_list[-i] <= vt_list[-1] - last_vi_margin) or (vt_list[-i] >= vt_list[-1] + last_vi_margin):
                return
        self.steady_state = True
        return

    def create_histogram(self):
        """
        Running ions in the system and collecting their location after the simulation.
        Creating a histogram from all collected locations
        """
        x_results = []
        number_of_cycles_per_ion = 10
        simulation_count = 0
        ion_subject = self.create_ion()
        if not self.fast_mode:
            print("\nRIMS simulation in progress...")
            create_trace_file(self, ion_subject)
        vT_list = []
        '''main simulation loop'''
        # while not self.steady_state:
        vt_list = []
        #     self.steady_state_matrix = []
        while simulation_count < self.number_of_simulations:
            ion_subject = self.create_ion()
            x_results.append(ion_subject.simulate_ion(number_of_cycles_per_ion))
            self.steady_state_matrix.append(ion_subject.velocity_list)
            # vt_list.append(ion_subject.steady_state_velocity)
            simulation_count += 1

            if not self.fast_mode:
                percentage_progress(simulation_count, self.number_of_simulations)
                write_to_trace_file(self, ion_subject)
        for j in range(POINTS):
            vtj_list = []
            for k in range(self.number_of_simulations):
                vtj_list.append(self.steady_state_matrix[k][j])
            vtj_array = np.array(vtj_list)
            vtj_av_speed = np.average(vtj_array)
            vt_list.append(vtj_av_speed)

        v_plot_list = []
        #rims.check_for_steady_state(self, vt_list)

        for v in range(len(vt_list)):
            if (v % 13 == 0):
                v_plot_sliced_array = np.array(vt_list[v - 13 : v])
                v_plot_list.append(np.average(v_plot_sliced_array))

        unique_id = create_unique_id()
        plt.figure(unique_id)
        x_axis = [cycle + 1 for cycle in range(len(v_plot_list))]

        plt.plot(x_axis, v_plot_list)
        plt.xlabel(r"Ratchet Cycle")
        plt.ylabel(r"Particle Velocity [cm/sec]")
        plt.suptitle("Average speed of ions over ratchet cycles")
        plt.savefig("Average speed of ions over ratchet cycles.jpeg")
        plt.close(unique_id)

        # '''calculation of particles velocity and current at steady state'''
        # vt_vector = np.array(vt_list)
        # vt_av = np.average(vt_vector)
        # vt_over_T = vt_av * ion_subject.interval * self.flash_frequency
        # vT_list.append(vt_av)
        # rims.check_for_steady_state(self, vT_list)
        vT_array = np.array(vt_list[-5:])
        vT_av = np.average(vT_array)
        number_of_cycles_per_ion = min(number_of_cycles_per_ion * 2, int(POINTS/14))
        simulation_count = 0
        print("number of cycles per ion is : " + str(number_of_cycles_per_ion))
        self.current = get_current(-vT_av, NE, SIGMA, ELECTRON_CHARGE)
        print("number of cycles per ion is : " + str(number_of_cycles_per_ion))

        if not self.fast_mode:
            print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
            create_log_file(self, ion_subject)
            print("Simulation log file created and saved.\n")

            '''plotting distribution of particles histogram'''
            plot_id = create_unique_id()
            plt.figure(plot_id)
            weights = np.ones_like(x_results) / float(len(x_results))
            x_results_um = [x * np.power(10, 4) for x in x_results]
            plt.hist(x_results_um, weights=weights, bins=RESOLUTION, label=r"X [$\mu $m]")
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


'''----------------------------------------------------------------------
                            CLASS IMPLEMENTATION
----------------------------------------------------------------------'''


def execution(dc_sample, f_sample, fast_mode):
    """
    :param dc_sample: duty cycle of positive part of the potential profile (only for fast mode)
    :param f_sample: frequency of the ratchet (only for fast mode)
    :param fast_mode: no logs or plots of execution will be saved, no prints to console
    creating a new simulation environment and launching it
    """
    if fast_mode:
        '''no logs or plots of execution will be saved, no prints to console'''
        ion_selection = debug_dict["ion_selection"]
        potential_profile = debug_dict["potential_profile"]
        flash_frequency = f_sample
        dc = dc_sample
        flash_mode = debug_dict["flash_mode"]
        output_selection = debug_dict["output_selection"]

    else:
        number_selection = ion_selection_panel()

        if number_selection == 6:
            '''Debug mode: preselected parameters for simulator functional testing'''
            ion_selection = debug_dict["ion_selection"]
            potential_profile = debug_dict["potential_profile"]
            flash_frequency = 600000
            dc = 0.6
            flash_mode = debug_dict["flash_mode"]
            output_selection = debug_dict["output_selection"]

        else:
            '''Extraction of all the data from the user'''
            ion_selection, potential_profile, flash_frequency, flash_mode, dc, output_selection = \
                extract_data_from_interface(number_selection)

    r = rims(ion_selection, potential_profile, flash_frequency, flash_mode, dc, fast_mode)
    r.get_electric_field()
    if output_selection == 1:
        r.create_histogram()
    elif output_selection == 2:
        r.create_video()
        generate_video_from_frames(r.path_for_output + 'frames', 'density over time.avi')

    if not fast_mode:
        print_log_file(r)
        rerun = execution_rerun_panel()
        if rerun:
            execution(dc_sample, f_sample, False)
    return r


def create_i_of_dc(frequency):
    """
    Runs 40 different dc samples to create I(dc) graph for constant frequency
    """
    dc_list = []
    current_list = []
    plot_uid = create_unique_id()
    runs = 40
    start_time = datetime.now()
    print("Creating I(duty cycle) graph for frequency="+str(frequency))
    for dc_sample in range(0, runs):
        percentage_progress(dc_sample, runs)
        rims_object = execution(dc_sample/runs, frequency, True)
        current_calculated = rims_object.current
        dc_list.append(dc_sample/runs)
        current_list.append(current_calculated)
    percentage_progress(1, 1)

    '''Plotting graph, I as a function of DC for constant frequency'''
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")
    plt.figure(plot_uid)
    plt.plot(dc_list, current_list, label="I(duty cycle)", color='#f5bc42')
    plt.suptitle('RIMS: current changing over DC at '+str(frequency)+'Hz', fontsize=12, fontweight='bold')
    plt.xlabel(r"DC")
    plt.ylabel(r"I [AMP]")

    '''Adding min & max markers'''
    max_current = max(current_list)
    max_i = current_list.index(max_current)
    plt.plot(dc_list[max_i], max_current, 'g^')
    plt.text(dc_list[max_i], max_current, str(max_current))

    min_current = min(current_list)
    min_i = current_list.index(min_current)
    plt.plot(dc_list[min_i], min_current, 'rv')
    plt.text(dc_list[min_i], min_current, str(min_current))

    '''Zero current line'''
    plt.axhline(color='r')
    file_name = 'i_dc ' + plot_uid + '.jpeg'
    plt.savefig(file_name)
    print('I_dc saved to main directory as ' + file_name)
    return


def create_heat_map():
    """
    Creates a heat-map of current as a function of DC anf frequency
    """
    resolution_f = 10           # number of frequencies tested
    resolution_dc = 10          # number of Dc's tested
    dc_list = [dc/resolution_dc for dc in range(resolution_dc, 0, -1)]
    f_list = [f * 1000 for f in range(100, 1100, 100)]

    plot_uid = create_unique_id()
    matrix = np.zeros(shape=(resolution_f, resolution_dc))

    start_time = datetime.now()
    for i_f in range(0, resolution_f):
        current_vector = np.zeros(resolution_dc)
        for i_dc in range(0, resolution_dc):
            percentage_progress(i_dc + i_f * resolution_dc, resolution_dc * resolution_f)
            dc_input = dc_list[i_dc]
            f_input = f_list[i_f]
            current_calculated = execution(dc_input, f_input, True).current
            current_vector[i_dc] = current_calculated

        matrix[resolution_f - i_f - 1] = current_vector
    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")

    '''Plotting heat-map'''
    fig, ax = plt.subplots()
    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    im, _ = heatmap(matrix.transpose(), dc_list, f_list_label, ax=ax, vmin=np.min(matrix),
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


'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''

headline_panel()
# create_i_of_dc(frequency=200000)
# create_heat_map()
execution(0, 0, False)
