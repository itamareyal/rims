from ion_simulation import *
from current_calc import *
from output import *
from defines import *
from interface import *
import PySimpleGUIQt as sg

""""
rims.py

---TOP MODULE---

Simulator of the movement of ions in a ratchet based system
for further information regarding this project, refer to README file.

This script hosts the Rims class and defining simulation parameters. 
In addition, it holds horizontal execution functions.

Eran Weil
Itamar Eyal
Dr. Gideon Segev
Energy devices lab, Tel-Aviv university
"""

'''----------------------------------------------------------------------
                         CLASS IMPLEMENTATION
----------------------------------------------------------------------'''


class Rims:
    def __init__(self, ion_subject, potential_profile, fast_mode):
        """
        Instance holding attributes of the whole simulation
        :param ion_subject: ion type to be simulated
        :param potential_profile: list holding potential profile shape [L,a1,a2]
        :param fast_mode: indication for no plotting and mid-calculation data savings
        """

        '''ratchet attributes'''
        self.ion = ion_subject[0]
        self.start_time = datetime.now()
        self.path_for_output = r'simulation outputs/' + get_time_stamp(self.start_time) + ' ' + self.ion + r'/'
        self.log = create_log_file(self)
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
        self.settings = load_settings(settings_filename)
        self.fast_mode = fast_mode
        self.steady_state = False
        self.resolution = RESOLUTION if self.potential_profile_mat.shape[1] < RESOLUTION\
            else self.potential_profile_mat.shape[1]
        self.electric_field_mat = self.get_electric_field()
        self.ions_lst = [Rims.create_ion(self) for _ in range(self.settings['PARTICLES_SIMULATED'])]
        self.gui = None

        '''result attributes'''
        self.x_results = []
        self.frames = np.zeros(shape=(load_one_setting(settings_filename,'MAX_CYCLES')+1, self.settings['PARTICLES_SIMULATED']+1))
        self.velocity = 0
        self.var = 0
        self.current = 0
        self.css = 0
        write_to_log(self, "RIMS simulation initialized")

    def create_ion(self):
        """
        Creates ion class instance to be simulated
        :return: ion instance
        """
        return Ion(self)

    def get_electric_field(self):
        """
        Derives the electric field from the potential, E(x,t) saves it as attribute
        """
        dx = self.L/self.resolution
        self.electric_field_mat = np.array([-np.gradient(v, dx) for v in self.potential_profile_mat])
        '''plot E & V'''
        write_to_log(self, "electric field generated")
        if not self.fast_mode:
            plot_potential_profile(self)
        return self.electric_field_mat

    def get_intervals(self):
        """
        Calculates the length of a time interval delta_t
        """
        '''Check for manual dt overwrite in settings'''
        if OVERWRITE_DELTA_T:
            write_to_log(self, "manual overwrite for delta t detected and is "+str(DELTA_T)+"[sec]")
            return DELTA_T
        '''Special dt for electrons'''
        if self.ion == "Electrons in Silicon":
            return self.flash_period / INTERVALS_FLASH_RATIO_ELECTRONS
        '''keeps delta t smaller than tau (diffusion time for L)'''
        critical_t = ((1 / INTERVALS_FLASH_RATIO) * self.L) ** 2 / (2 * self.diffusion)
        '''Further diminishes delta t if its bigger that T/INTERVALS_FLASH_RATIO'''
        while critical_t > self.flash_period / INTERVALS_FLASH_RATIO:
            critical_t /= INTERVALS_FLASH_RATIO
        write_to_log(self, "delta t calculated and is "+str(critical_t)+"[usec]")
        return critical_t

    def get_gamma(self):
        gamma = BOLTZMANN_CONSTANT * TEMPERATURE / self.diffusion
        write_to_log(self, "gamma factor calculated "+str(gamma)+"[eVsec/cm^2]")
        return gamma

    def check_for_steady_state(self, vt_list):
        """
        Checks whether the system has reached ss
        """
        num_discrepancies_allowed = 1
        '''Allow for ss only after MIN_MEASUREMENTS_FOR_SS cycles'''
        if len(vt_list) < MIN_MEASUREMENTS_FOR_SS:
            write_to_log(self, "cycle count still lower than MIN_MEASUREMENTS_FOR_SS="+str(MIN_MEASUREMENTS_FOR_SS))
            return
        margin = abs(vt_list[-1]) * self.settings['STEADY_STATE_PERCENT_MARGIN']
        mean = np.average(vt_list[-5:])
        discrepancies = 0
        for i in range(5):
            if (mean-margin) <= vt_list[i] <= (mean+margin):
                continue
            '''one discrepancy counted'''
            discrepancies += 1

        if discrepancies <= num_discrepancies_allowed:
            if not self.steady_state:
                '''Saves the cycle when ss was first reached'''
                self.css = self.cycles_count
                write_to_log(self, "steady state reached after "+str(self.css)+" cycles")
            else:
                write_to_log(self, "steady state maintained after "+str(self.css)+" cycles")
            self.steady_state = True
        else:
            write_to_log(self, "velocity is NOT within margin for steady state. discrepancies="+str(discrepancies))
        return

    def generate_ions_mat(self):
        ions_mat = []
        for thread in range(10):
            ions_lst = [Rims.create_ion(self) for _ in range(IONS_PER_THREAD)]
            ions_mat.append(ions_lst)
        return ions_mat

    def get_velocity_over_cycle(self):
        """
        Runs all ions for a single ratchet cycle and calculates averaged velocity over that cycle
        """
        cycle_v = []
        for x, ion_subject in enumerate(self.ions_lst):
            '''simulate all ions for 1 cycle'''
            ion_subject.simulate_ion()
            cycle_v.append(ion_subject.velocity)
            if load_one_setting(settings_filename,'ENABLE_VIDEO'):
                self.frames[self.cycles_count][x] = ion_subject.absolute_final_loc  # collect for video
        cycle_v = np.array(cycle_v)
        write_to_log(self, "cycle " + str(self.cycles_count) + " completed. " +
                                                               "velocity is " + str(np.average(cycle_v)) + "[cm/sec]")
        self.gui.Refresh() if self.gui else None
        return np.average(cycle_v)

    def run_rims(self):
        """
        Running ions in the system and collecting their location after the simulation.
        Creating a histogram from all collected locations
        """
        if not self.fast_mode:
            print("\nRIMS simulation in progress...")

        '''Holds average velocity of all ions per cycle. rims_v[i]=av velocity at cycle i'''
        rims_v = []

        '''Main simulation loop'''
        while self.cycles_count <= load_one_setting(settings_filename,'MAX_CYCLES'):
            #sg.one_line_progress_meter('RIMS simulation', self.cycles_count, MAX_CYCLES, 'key',
            #'Ion type: '+self.ion, orientation='h')

            percentage_progress(self.cycles_count, load_one_setting(settings_filename,'MAX_CYCLES'), self)
            rims_v.append(self.get_velocity_over_cycle())
            self.check_for_steady_state(rims_v)
            self.cycles_count += 1

        '''Collecting final locations for histogram'''
        for ion_subject in self.ions_lst:
            self.x_results.append(ion_subject.absolute_final_loc)

        write_to_log(self, "data collected for histogram")

        self.velocity = np.average(rims_v[1:])
        self.var = np.var(rims_v[1:])
        self.current = get_current(self.velocity, NE, SIGMA, ELECTRON_CHARGE)

        if not self.fast_mode:
            print("\nSimulation finished after " + str(datetime.now() - self.start_time) + "\n")
            create_summary_file(self)
            
            create_trace_file(self) if load_one_setting(settings_filename,'CREATE_TRACE') else None

            '''plotting distribution of particles histogram & average speed plot'''
            x_results = np.array(self.x_results)
            plot_distribution_over_x_histogram(self, x_results)
            plot_average_speed_of_ions(self, rims_v)
            print("All plots and logs saved to folder " + self.path_for_output + '\n')
        write_to_log(self, "simulation finished")
        return


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''

def thread_main(ions_lst, v_cycle):
    thread_v = []
    for ion_subject in ions_lst:
        ion_subject.simulate_ion()  # simulate for 1 cycle
        thread_v.append(ion_subject.velocity)
    thread_v = np.array(thread_v)
    v_cycle.append(np.average(thread_v))
    return

def set_gui(rims_object, ions_simulated, pb_iter, gui):
    gui.Refresh() if gui else None
    rims_object.gui = gui
    if gui:
        rims_object.pb = gui['progressbar']
        rims_object.pb_iter = pb_iter
        rims_object.pb_max = ions_simulated * load_one_setting(settings_filename,'MAX_CYCLES')

def execution(ion_selection_dict, potential_profile, gui=None):
    """
    creating a new simulation environment and launching it
    """
    headline_panel()
    simulations = []

    video_3d_mat = []

    pb_iter = 0

    '''Simulating for every ion specified'''
    for key, value in ion_selection_dict.items():
        ion_selection = (key, value)
        print('\n-------------------------------------------------------\n')
        print('Simulating '+ion_selection[0]+'; D='+str(ion_selection[1])+'[cm^2/sec]')
        r = Rims(ion_selection, potential_profile, False)
        set_gui(r, len(ion_selection_dict.items()), pb_iter, gui)
        r.run_rims()
        simulations.append(r)
        pb_iter = r.pb_iter if gui else None

        '''Video'''
        video_3d_mat.append(r.frames) if load_one_setting(settings_filename,'ENABLE_VIDEO') else None

        print_summary_file(r)

    plot_multiple_ions_hist(simulations)
    '''Video'''
    if load_one_setting(settings_filename,'ENABLE_VIDEO'):
        create_video_of_histograms(video_3d_mat, ion_selection_dict, simulations[0])

    simulation_over_panel()
    if __name__ == '__main__':
        if execution_rerun_panel():
            execution(ion_selection_dict, potential_profile)

    return

def plot_multiple_ions_hist(simulations):
    plot_id = create_unique_id()
    plt.figure(plot_id)
    for r in simulations:
        '''Adding distribution data for the complete histogram'''
        x = r.x_results
        weights = np.ones_like(x) / float(len(x))
        x_um = [x * np.power(10, 4) for x in x]
        ion_selection = (r.ion, r.diffusion)
        if ion_selection[0][:6] == 'manual':
            label = 'D='+str(ion_selection[1])
        else:
            label = ion_selection[0]
        label += '; '+"{:.3f}".format(r.velocity)+'[cm/sec]'
        plt.hist(x_um, weights=weights, bins=r.resolution, label=label)
        write_to_log(r, "Data added to multiple ions histogram")

    '''Plotting combined histogram of ions simulated'''
    plt.ylabel('Density')
    plt.xlabel(r'X [$\mu $m]')
    plt.legend()
    plt.title(r"RIMS: Histogram of distribution x axis: $\rho $(x)", fontsize=14, fontweight='bold')
    time_stamp = get_time_stamp(datetime.now())
    file_name = str(time_stamp) + ' Distribution histogram'

    for r in simulations:
        if not os.path.exists(r.path_for_output):
            os.makedirs(r.path_for_output)
        plt.savefig(r.path_for_output + file_name)
        write_to_log(r, "Multiple ions histogram saved")

    folder = r'simulation outputs/00_Multiple ions histogram'
    if not os.path.exists(folder):
        os.makedirs(folder)
    plt.savefig(folder+r'/'+file_name)
    plt.close(plot_id)
    print('Histogram of all ions simulated saved in '+folder+' as '+file_name+'.jpeg')



'''----------------------------------------------------------------------
                               EXECUTION
----------------------------------------------------------------------'''

if __name__ == '__main__':

    '''Extraction of data from interface'''
    ion_selection_dict = ion_selection_panel()
    potential_profile = extract_data_from_interface()
    '''Running simulation'''
    execution(ion_selection_dict, potential_profile)
