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
    def __init__(self, ion_subject, potential_profile, fast_mode):
        """
        Instance holding attributes of the whole simulation
        :param ion_subject: ion type to be simulated
        :param potential_profile: list holding potential profile shape [L,a1,a2]
        :param flash_frequency: flashing rate of the ratchet
        :param fast_mode: indication for no plotting and mid-calculation data savings
        """
        '''ratchet attributes'''
        self.ion = ion_subject[0]
        self.current = 0
        self.interval = -1
        self.electric_field = 0
        self.potential_profile = 0
        self.potential_profile_list = potential_profile
        self.diffusion = ion_subject[1]
        self.L = self.potential_profile_list[0]                         # profile length x axis [um]
        self.x_space_vec = self.potential_profile_list[1]
        self.time_vec = self.potential_profile_list[2]
        self.potential_profile_mat = self.potential_profile_list[3]
        self.electric_field_mat = np.zeros(shape=(1, 1))                # initialized to zeros
        self.flash_period = self.time_vec[-1]
        self.flash_frequency = 1/self.time_vec[-1]

        '''simulation attributes'''
        self.number_of_simulations = NUMBER_OF_SIMULATIONS
        self.start_time = datetime.now()
        self.path_for_output = r'RIMS output plots/' + get_time_stamp(self.start_time) + ' ' + self.ion + r'/'
        self.fast_mode = fast_mode
        self.steady_state_matrix = []
        self.steady_state = False
        self.resolution = RESOLUTION if self.potential_profile_mat.shape[1] < RESOLUTION\
            else self.potential_profile_mat.shape[1]

        '''result attributes'''
        self.velocity = 0
        self.current = 0

    def create_ion(self):
        """
        Creates ion class instance to be simulated
        :return: ion instance
        """
        i = ion(self.diffusion, self.L, self.time_vec, self.electric_field_mat, self.path_for_output)
        i.get_intervals()
        i.get_gamma()
        return i

    def get_electric_field(self):
        """
        Derives the electric field from the potential, E(x,t) saves it as attribute
        """

        self.electric_field_mat = np.array([-np.gradient(v, self.L/self.resolution) for v in self.potential_profile_mat])

        '''plot E & V'''
        if not self.fast_mode:
            plot_potential_profile(self)
        return


    def check_for_steady_state(self, vt_list):
        num_discrepancies_allowed = 1
        if len(vt_list) < MIN_NUM_SPEEDS_FOR_AVG:
            return
        last_vi_margin = vt_list[-1] / 10
        for i in range(2, 6, 1):
            if (vt_list[-i] <= vt_list[-1] - last_vi_margin) or (vt_list[-i] >= vt_list[-1] + last_vi_margin):
                if num_discrepancies_allowed > 0:
                    num_discrepancies_allowed -= 1
                else:
                    return

        self.steady_state = True
        return

    def create_histogram(self):
        """
        Running ions in the system and collecting their location after the simulation.
        Creating a histogram from all collected locations
        """
        number_of_cycles_per_ion = 1
        ion_subject = self.create_ion()
        v_plot_list = []
        x_results = []

        if not self.fast_mode:
            print("\nRIMS simulation in progress...")
            create_trace_file(self, ion_subject)

        '''main simulation loop'''
        while not self.steady_state:
            x_results = []
            vt_list = []
            self.steady_state_matrix = []
            simulation_count = 0
            while simulation_count < self.number_of_simulations:
                ion_subject = self.create_ion()
                x_results.append(ion_subject.simulate_ion(number_of_cycles_per_ion))
                self.steady_state_matrix.append(ion_subject.velocity_list)
                simulation_count += 1

                if not self.fast_mode:
                    percentage_progress(simulation_count, self.number_of_simulations)
                    write_to_trace_file(self, ion_subject)

            for j in range(len(self.steady_state_matrix[0])):
                vtj_list = []
                for k in range(self.number_of_simulations):
                    vtj_list.append(self.steady_state_matrix[k][j])
                vtj_array = np.array(vtj_list)
                vtj_av_speed = np.average(vtj_array)
                vt_list.append(vtj_av_speed)

            v_plot_list = []
            num_of_intervals_per_cycle =int(self.flash_period / ion_subject.interval)
            for v in range(len(vt_list)):
                if (v % num_of_intervals_per_cycle == 0) and (v != 0):
                    v_plot_sliced_array = np.array(vt_list[v - num_of_intervals_per_cycle : v])
                    v_plot_list.append(np.average(v_plot_sliced_array))

            rims.check_for_steady_state(self, v_plot_list)
            self.steady_state = True
            #number_of_cycles_per_ion = number_of_cycles_per_ion * 2

        '''calculation of particles velocity and current at steady state'''
        if len(v_plot_list) >= MIN_NUM_SPEEDS_FOR_AVG:
            vT_av_array = np.array(v_plot_list[-MIN_NUM_SPEEDS_FOR_AVG:])
            vT_av = np.average(vT_av_array)
            self.velocity = vT_av
            if not self.fast_mode:
                plot_average_speed_of_ions(self, v_plot_list)
            self.current = get_current(vT_av, NE, SIGMA, ELECTRON_CHARGE)

        # number_of_cycles_per_ion = min(number_of_cycles_per_ion * 2, int(POINTS/14))
        # print("number of cycles per ion is : " + str(number_of_cycles_per_ion))

        if not self.fast_mode:
            print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
            create_log_file(self, ion_subject)
            print("Simulation log file created and saved.\n")

            '''plotting distribution of particles histogram'''
            x_results = np.array(x_results)
            plot_distribution_over_x_histogram(self, x_results)
        return x_results

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
            plt.hist(x_results, weights=weights, bins=self.resolution, label=r"X [$\mu $m]", range=(-5, 5))
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
        '''Fast mode: no logs or plots of execution will be saved, no prints to console'''
        ion_selection = debug_dict["ion_selection"]
        potential_profile = debug_dict["potential_profile"]
        potential_profile[2][0] = float(dc_sample / f_sample)
        potential_profile[2][1] = 1 / f_sample
        r = rims(ion_selection, potential_profile, fast_mode)
        r.get_electric_field()
        r.create_histogram()

    else:
        ion_selection_dict = ion_selection_panel()

        if ion_selection_dict == 6:
            '''Debug mode: preselected parameters for simulator functional testing'''
            ion_selection = (debug_dict["ion_selection"])
            potential_profile = debug_dict["potential_profile"]
            r = rims(ion_selection, potential_profile, fast_mode)
            r.get_electric_field()
            r.create_histogram()
            print('i='+str(r.current))

        else:
            '''Extraction of all the data from the user'''
            potential_profile = extract_data_from_interface()

            '''Simulating for every ion specified'''
            plot_id = create_unique_id()
            plt.figure(plot_id)
            for key, value in ion_selection_dict.items():
                ion_selection = (key, value)
                print('\n-------------------------------------------------------\n')
                print('Simulating '+ion_selection[0]+'; Deff='+str(ion_selection[1])+'[cm^2/sec]')
                r = rims(ion_selection, potential_profile, fast_mode)
                r.get_electric_field()
                x_results = r.create_histogram()
                weights = np.ones_like(x_results) / float(len(x_results))
                x_results_um = [x * np.power(10, 4) for x in x_results]
                if ion_selection[0][:7] == 'manual':
                    label = str(ion_selection[1])
                else:
                    label = ion_selection[0]
                plt.hist(x_results_um, weights=weights, bins=r.resolution, label=label)
                print_log_file(r)

            '''Plotting combined histogram of ions simulated'''
            plt.ylabel('Density')
            plt.xlabel(r'X [$\mu $m]')
            plt.legend()
            plt.title(r"RIMS: Histogram of distribution x axis: $\rho $(x)", fontsize=14, fontweight='bold')
            time_stamp = get_time_stamp(datetime.now())
            file_name = str(time_stamp) + ' Distribution histogram'
            plt.savefig(file_name)
            plt.close(plot_id)
            print('Histogram saved as '+file_name+'.jpeg')

    if not fast_mode:
        rerun = execution_rerun_panel()
        if rerun:
            execution(dc_sample, f_sample, False)
    return r


def create_i_of_dc_comparison(frequencies, compare):
    """
    Runs different dc samples to create I(dc) graph for constant frequency and compare with Kedem's exp
    :param frequencies: list of f to be simulated over every dc
    :param compare: bool, also plot analytical calculation by Kedem
    """
    plot_uid = create_unique_id()
    runs = 20
    dc_list = [dc/runs for dc in range(0,runs+1)]
    start_time = datetime.now()
    print("Creating I(duty cycle) graph")
    plt.figure(plot_uid)
    i_f = 0
    for frequency in frequencies:
        current_list = []
        k_currents = []
        for dc_sample in range(0, runs+1):
            dc = dc_sample/runs
            percentage_progress(dc_sample + i_f*len(frequencies), runs * len(frequencies))
            rims_object = execution(dc, frequency, True)
            current_calculated = rims_object.current
            current_list.append(current_calculated)

            if compare:
                '''collect for Kedem'''
                k_v = get_velocity(float(1/frequency), 0.8 * pow(10, -4),
                                   2.5 * pow(10, -4),
                                   0.25*ELECTRON_CHARGE, 0.05*ELECTRON_CHARGE, -0.5,
                                   TEMPERATURE, 1 - float(dc))

                k_i = get_current(k_v, NE, SIGMA, ELECTRON_CHARGE)
                k_currents.append(k_i/10)
        i_f += 1
        plt.plot(dc_list, current_list, label='RIMS: '+str(frequency / 1000) + "KHz")
        if compare:
            plt.plot(dc_list, k_currents, label='Kedem: '+str(frequency / 1000) + "KHz")
    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")
    plt.suptitle('RIMS: current changing over DC', fontsize=12, fontweight='bold')
    plt.xlabel(r"DC")
    plt.ylabel(r"I [AMP]")
    plt.legend(loc='upper left')

    '''Zero current line'''
    plt.axhline(color='r')
    file_name = 'i_dc changing f' + plot_uid + '.jpeg'
    plt.savefig(file_name)
    print('I_dc saved to main directory as ' + file_name)
    plt.close(plot_uid)
    return

def create_i_of_f_comparison(dc, compare):
    """
    Runs different frequency samples to create I(f) graph for constant dc and compare with Kedem's exp
    """
    f_list = [f * 1000 for f in range(10, 1100, 100)]
    current_list = []
    k_currents = []
    plot_uid = create_unique_id()
    start_time = datetime.now()
    print("Creating I(f) graph for dc="+str(dc))
    percentage_progress(0,1)
    for f_sample in f_list:
        rims_object = execution(dc, f_sample, True)
        current_calculated = rims_object.current
        current_list.append(current_calculated)

        if compare:
            '''collect for Keden'''
            k_v = get_velocity(1/f_sample, 0.8 *pow(10,-4),
                               1.3 * pow(10, -5),
                               0.25*ELECTRON_CHARGE, 0.05*ELECTRON_CHARGE, -0.5, TEMPERATURE,1-dc)
            k_i = get_current(k_v, NE, SIGMA, -ELECTRON_CHARGE)
            k_currents.append(k_i)
        percentage_progress(f_sample, max(f_list))

    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")

    '''Plotting graph, I as a function of DC for constant frequency'''
    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    plt.figure(plot_uid)
    plt.plot(f_list_label, current_list, label="RIMS", color=PURPLE)
    if compare:
        plt.plot(f_list_label, k_currents, label='Kedem exp', color=YELLOW)
    plt.suptitle('RIMS: current changing over f at dc='+str(dc), fontsize=12, fontweight='bold')
    plt.xlabel(r"f [KHz]")
    plt.ylabel(r"I [AMP]")
    plt.legend(loc='upper right')

    '''Adding min & max markers'''
    max_current = max(current_list)
    max_i = current_list.index(max_current)
    plt.plot(f_list_label[max_i], max_current, 'g^')
    plt.text(f_list_label[max_i], max_current, str(max_current))

    min_current = min(current_list)
    min_i = current_list.index(min_current)
    plt.plot(f_list_label[min_i], min_current, 'rv')
    plt.text(f_list_label[min_i], min_current, str(min_current))

    '''Zero current line'''
    plt.axhline(color='r')
    file_name = 'i_f ' + plot_uid + '.jpeg'
    plt.savefig(file_name)
    print('I_f saved to main directory as ' + file_name)
    return

def create_heat_map():
    """
    Creates a heat-map of current as a function of DC anf frequency
    """
    resolution_f = 10           # number of frequencies tested
    resolution_dc = 10          # number of Dc's tested
    dc_list = [dc/resolution_dc for dc in range(resolution_dc, 0, -1)]
    f_list = [f * 1000 for f in range(10, 1010, 100)]

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

        matrix[:, i_f] = current_vector
    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")

    '''Plotting heat-map'''
    fig, ax = plt.subplots()
    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    bar_maximal_value = max(np.abs(np.min(matrix)), np.abs(np.max(matrix)))
    im, _ = heatmap(matrix, dc_list, f_list_label, ax=ax, vmin=-bar_maximal_value, vmax=bar_maximal_value,
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

# create_i_of_dc_comparison([800000,600000], False)
# create_i_of_f_comparison(0.6, True)
# create_heat_map()
execution(0, 0, False)
