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
import threading

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
        self.electric_field = 0
        self.potential_profile = 0
        self.diffusion = ion_subject[1]
        self.gamma = self.get_gamma()
        self.L = potential_profile[0]
        self.x_space_vec = potential_profile[1]
        self.time_vec = potential_profile[2]
        self.potential_profile_mat = potential_profile[3]
        self.flash_period = self.time_vec[-1]
        self.flash_frequency = 1/self.time_vec[-1]
        self.interval = self.get_intervals()
        self.intervals_in_period = int(self.flash_period/self.interval)
        self.cycles_count = 0

        '''simulation attributes'''
        self.number_of_simulations = NUMBER_OF_SIMULATIONS
        self.start_time = datetime.now()
        self.path_for_output = r'RIMS output plots/' + get_time_stamp(self.start_time) + ' ' + self.ion + r'/'
        self.fast_mode = fast_mode
        self.steady_state_matrix = []
        self.steady_state = False
        self.resolution = RESOLUTION if self.potential_profile_mat.shape[1] < RESOLUTION\
            else self.potential_profile_mat.shape[1]
        # self.ions_mat = self.generate_ions_mat()
        self.electric_field_mat = self.get_electric_field()
        self.ions_lst = [rims.create_ion(self) for i in range(NUMBER_OF_SIMULATIONS)]
        '''result attributes'''
        self.x_results = []
        self.velocity = 0
        self.current = 0

    def create_ion(self):
        """
        Creates ion class instance to be simulated
        :return: ion instance
        """
        return ion(self)

    def get_electric_field(self):
        """
        Derives the electric field from the potential, E(x,t) saves it as attribute
        """
        self.electric_field_mat = np.array([-np.gradient(v, self.L/self.resolution) for v in self.potential_profile_mat])
        '''plot E & V'''
        if not self.fast_mode:
            plot_potential_profile(self)
        return self.electric_field_mat

    def get_intervals(self):
        if self.ion == "Electrons in Silicon":
            return self.flash_period / 50
        critical_t = ((1 / INTERVALS_FLASH_RATIO) * self.L) ** 2 / (2 * self.diffusion)
        while critical_t > self.flash_period / INTERVALS_FLASH_RATIO:
            critical_t /= INTERVALS_FLASH_RATIO
        return critical_t

    def get_gamma(self):
        return BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion

    def check_for_steady_state(self, vt_list):
        num_discrepancies_allowed = 1
        if len(vt_list) < MIN_NUM_SPEEDS_FOR_AVG:
            return
        margin = abs(vt_list[-1]) * 0.1
        mean = np.average(vt_list[-5:])
        discrepancies = 0
        for i in range(5):
            if (mean-margin) <= vt_list[i] <= (mean+margin):
                continue
            discrepancies += 1

        if discrepancies <= num_discrepancies_allowed:
            self.steady_state = True
        return

    def generate_ions_mat(self):
        ions_mat = []
        for thread in range(NUMBER_OF_THREADS):
            ions_lst = [rims.create_ion(self) for i in range(IONS_PER_THREAD)]
            ions_mat.append(ions_lst)
        return ions_mat

    def get_velocity_over_cycle(self):
        # threads = []
        cycle_v = []
        # thread_v = []
        # while cycle_index < num_of_cycles:
        for ion_subject in self.ions_lst:
            ion_subject.simulate_ion()  # simulate for 1 cycle
            cycle_v.append(ion_subject.velocity)

        return np.average(np.array(cycle_v))

    def run_rims(self):
        """
        Running ions in the system and collecting their location after the simulation.
        Creating a histogram from all collected locations
        """
        if not self.fast_mode:
            print("\nRIMS simulation in progress...")
            create_trace_file(self)

        rims_v = []

        '''Main simulation loop'''
        while (not self.steady_state) and (self.cycles_count < MAX_CYCLES):
            if not self.fast_mode:
                percentage_progress(self.cycles_count, MAX_CYCLES)
            rims_v.append(self.get_velocity_over_cycle())
            self.check_for_steady_state(rims_v)
            self.cycles_count += 1

        if not self.fast_mode:
            percentage_progress(1, 1)
            print('\n')
            if self.steady_state:
                print('Steady state reached after '+str(self.cycles_count)+' ratchet cycles')
            else:
                print('Steady state NOT reached after '+str(self.cycles_count)+' ratchet cycles')
                print('Maximal number of cycles can be edited in defines.py module')

        '''Collecting final locations for histogram'''
        for ion_subject in self.ions_lst:
            self.x_results.append(ion_subject.absolute_final_loc)
            if not self.fast_mode:
                write_to_trace_file(self, ion_subject)

        self.velocity = np.average(rims_v[1:])
        self.current = get_current(self.velocity, NE, SIGMA, ELECTRON_CHARGE)

        if not self.fast_mode:
            print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
            create_log_file(self)
            print("Simulation log file saved. plotting graphs...")

            '''plotting distribution of particles histogram & average speed plot'''
            x_results = np.array(self.x_results)
            plot_distribution_over_x_histogram(self, x_results)
            plot_average_speed_of_ions(self, rims_v)
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

            create_trace_file(self)

            while simulation_count < self.number_of_simulations:
                ion_subject = self.create_ion()
                ion_subject.points = frame
                for cycle in MAX_CYCLES:
                    ion_subject.simulate_ion()
                x_results.append(ion_subject.absolute_final_loc)
                simulation_count += 1
                percentage_progress(simulation_count * 100, NUMBER_OF_SIMULATIONS)
                write_to_trace_file(self, ion_subject)

            plt.figure(frame)
            weights = np.ones_like(x_results) / float(len(x_results))
            plt.hist(x_results, weights=weights, bins=self.resolution, label=r"X [$\mu $m]", range=(-5, 5))
            plt.ylim(0, 0.05)
            plt.ylabel('Density')
            plt.xlabel(r'X [$\mu $m]')
            plt.title(r"RIMS: Histogram of distribution x axis: $\rho $(x,t)", fontsize=12, fontweight='bold')
            plt.suptitle('t = ' + str(ion_subject.interval * ion_subject.intervals_count)[0:8] + r' [$\mu $sec]',
                         fontsize=10)
            save_plots(self, 'frames\\frame_' + str(frame), frame)
            plt.close(frame)
        print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
        create_log_file(self)
        print("Simulation log file created and saved.\n")


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''

def thread_main(ions_lst, v_cycle):
    thread_v = []
    # while cycle_index < num_of_cycles:
    for ion_subject in ions_lst:
        ion_subject.simulate_ion() # simulate for 1 cycle
        thread_v.append(ion_subject.velocity)
    thread_v = np.array(thread_v)
    v_cycle.append(np.average(thread_v))

def debug_execution():
    '''Debug mode: preselected parameters for simulator functional testing'''
    ion_selection = (debug_dict["ion_selection"])
    potential_profile = debug_dict["potential_profile"]
    r = rims(ion_selection, potential_profile, False)
    r.run_rims()
    print('i=' + str(r.current))

def execution():
    """
    creating a new simulation environment and launching it
    """

    '''Extraction of data from interface'''
    ion_selection_dict = ion_selection_panel()
    potential_profile = extract_data_from_interface()

    plot_id = create_unique_id()
    plt.figure(plot_id)

    '''Simulating for every ion specified'''
    for key, value in ion_selection_dict.items():
        ion_selection = (key, value)
        print('\n-------------------------------------------------------\n')
        print('Simulating '+ion_selection[0]+'; D='+str(ion_selection[1])+'[cm^2/sec]')
        r = rims(ion_selection, potential_profile, False)
        r.run_rims()

        '''Adding distribution data for the complete histogram'''
        x = r.x_results
        weights = np.ones_like(x) / float(len(x))
        x_um = [x * np.power(10, 4) for x in x]
        if ion_selection[0][:6] == 'manual':
            label = 'D='+str(ion_selection[1])
        else:
            label = ion_selection[0]
        label += '; '+"{:.3f}".format(r.velocity)+'[cm/sec]'
        plt.hist(x_um, weights=weights, bins=r.resolution, label=label)
        #print_log_file(r)

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

    rerun = execution_rerun_panel()
    if rerun:
        execution()
    return


def create_i_of_dc_comparison(frequencies, compare):
    """
    Runs different dc samples to create I(dc) graph for constant frequency and compare with Kedem's exp
    :param frequencies: list of f to be simulated over every dc
    :param compare: bool, also plot analytical calculation by Kedem
    """
    plot_uid = create_unique_id()
    runs = 15
    dc_list = [dc/runs for dc in range(0, runs+1)]
    start_time = datetime.now()
    print("Generating Velocity(duty cycle) graph...")
    plt.figure(plot_uid)
    i_f = 0
    for frequency in frequencies:
        current_list = []
        k_currents = []
        for i_dc in range(0, runs+1):
            dc = i_dc/runs
            percentage_progress(i_dc + i_f*len(frequencies), runs * len(frequencies))
            # ion_selection = debug_dict["ion_selection"]
            # potential_profile = debug_dict["potential_profile"]
            # potential_profile[2][0] = float(dc / frequency)
            # potential_profile[2][1] = 1 / frequency
            # r = rims(ion_selection, potential_profile, True)
            # r.run_rims()
            # current_list.append(r.current)

            if compare:
                '''collect for Kedem'''
                k_v = get_velocity(float(1/frequency), 0.8 * pow(10, -4),
                                   2.5 * pow(10, -4),
                                   0.25*ELECTRON_CHARGE, 0.05*ELECTRON_CHARGE, -1,
                                   TEMPERATURE, 1 - float(dc))

                k_i = get_current(k_v, NE, SIGMA, ELECTRON_CHARGE)
                k_currents.append(k_i*2/1000)
        i_f += 1
        # plt.plot(dc_list, current_list, label='RIMS: '+str(int(frequency / 1000)) + "KHz")
        if compare:
            plt.plot(dc_list, k_currents, label='Kedem: '+str(int(frequency / 1000)) + "KHz")

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
# debug_execution()
# execution()

create_i_of_dc_comparison([1000000], True)
