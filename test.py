from rims import *
"""
test.py

Test library for rims
"""

'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''

def debug_execution():
    """
    Execution in pre-defined parameters for system check
    """
    '''Debug mode: preselected parameters for simulator functional testing'''
    ion_selection = d_ion_selection
    potential_profile = d_potential_profile
    r = Rims(ion_selection, potential_profile, False)
    r.run_rims()
    print('v=' + str(r.velocity))

def test_delta_t(flash_ratio_lst):
    plot_uid = create_unique_id()
    runs = 5
    # dc_list = [dc/runs for dc in range(0, runs+1)]
    dc_list = [0.8]
    start_time = datetime.now()
    print("Generating Velocity(duty cycle) graph...")
    # fig, axs = plt.subplots(len(flash_ratio_lst), sharex='all', sharey='all')
    # fig.suptitle('Vertically stacked subplots')
    ion_selection = ("Lead Pb+2", diffusion_coefficient_dict["Lead Pb+2"])
    potential_profile_file = "Segev.csv"
    L, t_vec, potential_mat = load_data_from_csv("potential profile sources/"+potential_profile_file)
    x = np.linspace(0, L, num=potential_mat.shape[1])
    potential_profile = [L, x, t_vec, potential_mat]
    velocity_av = []
    var_lst = []
    simulation_mat = []

    for i_f, ratio in enumerate(flash_ratio_lst):
        velocity_lst = []
        simulation_lst = []
        local_var = []
        for i_run, run in enumerate(range(runs)):
            percentage_progress(i_run + i_f * runs, runs * len(flash_ratio_lst))

            r = Rims(ion_selection, potential_profile, True)
            r.interval = r.flash_period / ratio
            r.intervals_in_period = ratio
            r.ions_lst = [Rims.create_ion(r) for i in range(PARTICLES_SIMULATED)]
            r.run_rims()
            simulation_lst.append(r)
            velocity_lst.append(r.velocity)
            local_var.append(r.var)
        var_lst.append(np.average(local_var))
        velocity_av.append(np.average(velocity_lst))
        simulation_mat.append(simulation_lst)

    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")
    folder = 'delta_t_tests'
    if not os.path.exists(folder):
        os.makedirs(folder)

    '''documenting graph'''
    fig, axs = plt.subplots()
    axs.plot(flash_ratio_lst, var_lst)
    plt.suptitle("RIMS: velocity and variance as a function of T/dt", fontsize=12, fontweight='bold')
    plt.xlabel("T/dt")
    plt.ylabel("Variance")
    file_name = plot_uid + ' variance as a function of ratio'
    plt.savefig(folder + r'/' + file_name + '.jpeg')
    plt.close(plot_uid)

    with open(folder + r'/' + str(plot_uid)+" delta_t test log.txt", "a") as f:
        f.write("delta_t test log\n\n")
        f.write("\ttime created: " + str(plot_uid) + "\n")
        f.write("\ttest duration: " + str(datetime.now() - start_time) + "\n")
        f.write("\tion tested, Deff: " + str(ion_selection) + "\n")
        f.write("\tPotential profile file name: " + potential_profile_file + "\n")
        f.write("\n-------------RESULTS-------------\n")
        for i, sim in enumerate(simulation_mat):
            f.write("\tratio = "+str(flash_ratio_lst[i])+'\n')
            f.write("\tdelta_t = "+str(simulation_mat[i][0].flash_period / flash_ratio_lst[i])+'\n')
            f.write("\taverage var = "+str(var_lst[i])+'\n')
            f.write("\taverage speed = " + str(np.average([sim.velocity for sim in simulation_mat[i]])) + '[cm/sec]\n\n')
    f.close()
    print('Graph and log saved to '+folder)

    return

def create_i_of_dc_comparison(frequencies, compare):
    """
    Runs different dc samples to create I(dc) graph for constant frequency and compare with Kedem's exp
    :param frequencies: list of f to be simulated over every dc
    :param compare: bool, also plot analytical calculation by Kedem
    """
    plot_uid = create_unique_id()
    runs = 30
    dc_list = [dc/runs for dc in range(0, runs+1)]
    start_time = datetime.now()
    print("Generating Velocity(duty cycle) graph...")
    plt.figure(plot_uid)
    for i_f, frequency in enumerate(frequencies):
        velocity_list = []
        k_currents = []
        for i_dc, dc in enumerate(dc_list):
            percentage_progress(i_dc + i_f*len(frequencies), runs * len(frequencies))
            ion_selection = d_ion_selection
            potential_profile = d_potential_profile
            potential_profile[2][0] = float(dc / frequency)
            potential_profile[2][1] = 1 / frequency
            r = Rims(ion_selection, potential_profile, True)
            r.run_rims()
            velocity_list.append(r.velocity)

            if compare:
                '''collect for Kedem'''
                k_v = get_velocity(period=float(1/frequency),
                                   L=d_L,
                                   diffusion=d_diffusion,
                                   a1=d_a1*ELECTRON_CHARGE, a2=d_a2*ELECTRON_CHARGE, alpha=d_alpha,
                                   temperature=TEMPERATURE, dc=1-float(dc))
                k_currents.append(k_v/1000)
        plt.plot(dc_list, velocity_list, label='RIMS: '+str(int(frequency / 1000)) + "KHz")
        if compare:
            plt.plot(dc_list, k_currents, label='Kedem: '+str(int(frequency / 1000)) + "KHz")

    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")
    plt.suptitle('RIMS: current changing over DC', fontsize=12, fontweight='bold')
    plt.xlabel(r"DC")
    plt.ylabel(r"Particle velocity [cm/sec]")
    plt.legend(loc='upper left')
    plt.axhline(color='r')

    '''documenting graph'''
    folder = 'Velocity over duty cycle graphs'
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_name = plot_uid + ' Velocity over duty cycle'
    plt.savefig(folder + r'/' + file_name + '.jpeg')
    print('Graph saved to '+folder+' as ' + file_name)
    plt.close(plot_uid)
    return

def create_i_of_f_comparison(dc, compare):
    """
    Runs different frequency samples to create V(f) graph for constant dc and compare with Kedem's exp
    """
    f_list = [f * 1000 for f in range(10, 1100, 100)]
    velocity_list = []
    k_currents = []
    plot_uid = create_unique_id()
    start_time = datetime.now()
    print("Creating I(f) graph for dc="+str(dc))

    '''Iterating over different frequencies'''
    for i_f, frequency in enumerate(f_list):
        percentage_progress(i_f, len(f_list))
        ion_selection = d_ion_selection
        potential_profile = d_potential_profile
        potential_profile[2][0] = float(dc / frequency)
        potential_profile[2][1] = 1 / frequency
        r = Rims(ion_selection, potential_profile, True)
        r.run_rims()
        velocity_calculated = r.velocity
        velocity_list.append(velocity_calculated)

        if compare:
            '''collect for Kedem'''
            k_v = get_velocity(period=float(1 / frequency),
                               L=d_L,
                               diffusion=d_diffusion,
                               a1=d_a1 * ELECTRON_CHARGE, a2=d_a2 * ELECTRON_CHARGE, alpha=d_alpha,
                               temperature=TEMPERATURE, dc=1 - float(dc))
            k_currents.append(k_v)

    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")

    '''Plotting graph, I as a function of DC for constant frequency'''
    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    plt.figure(plot_uid)
    plt.plot(f_list_label, velocity_list, label="RIMS", color=PURPLE)
    if compare:
        plt.plot(f_list_label, k_currents, label='Kedem', color=YELLOW)
    plt.suptitle('RIMS: current changing over f at dc='+str(dc), fontsize=12, fontweight='bold')
    plt.xlabel(r"f [KHz]")
    plt.ylabel(r"Particle velocity [cm/sec]")
    plt.legend(loc='upper right')

    '''Adding min & max markers'''
    max_current = max(velocity_list)
    max_i = velocity_list.index(max_current)
    plt.plot(f_list_label[max_i], max_current, 'g^')
    plt.text(f_list_label[max_i], max_current, str(max_current))

    min_current = min(velocity_list)
    min_i = velocity_list.index(min_current)
    plt.plot(f_list_label[min_i], min_current, 'rv')
    plt.text(f_list_label[min_i], min_current, str(min_current))
    plt.axhline(color='r')

    '''documenting graph'''
    folder = 'Velocity over frequency graphs'
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_name = plot_uid + ' Velocity over frequency'
    plt.savefig(folder + r'/' + file_name + '.jpeg')
    print('Graph saved to '+folder+' as ' + file_name)
    plt.close(plot_uid)
    return

def create_heat_map():
    """
    Creates a heat-map of current as a function of DC anf frequency
    """
    resolution_f = 10           # number of frequencies tested
    resolution_dc = 10          # number of Dc's tested
    dc_list = [dc/resolution_dc for dc in range(resolution_dc, 0, -1)]
    f_list = [f * 1000 for f in range(100, 1000, 100)]
    plot_uid = create_unique_id()
    matrix = np.zeros(shape=(resolution_f, resolution_dc))
    start_time = datetime.now()
    '''Iterating over frequencies and duty cycles'''
    for i_f, f_input in enumerate(f_list):
        velocity_vector = np.zeros(resolution_dc)
        for i_dc, dc_input in enumerate(dc_list):
            percentage_progress(i_dc + i_f * resolution_dc, resolution_dc * resolution_f)
            dc_input = dc_list[i_dc]
            f_input = f_list[i_f]
            ion_selection = d_ion_selection
            potential_profile = d_potential_profile
            potential_profile[2][0] = float(dc_input / f_input)
            potential_profile[2][1] = 1 / f_input
            r = Rims(ion_selection, potential_profile, True)
            r.run_rims()
            velocity_calculated = r.velocity
            velocity_vector[i_dc] = velocity_calculated

        matrix[:, i_f] = velocity_vector
    percentage_progress(1, 1)
    print("\nSimulation finished after " + str(datetime.now() - start_time) + "\n")

    '''Plotting heat-map'''
    fig, ax = plt.subplots()
    f_list_label = [str(int(f / 1000)) + 'k' for f in f_list]
    bar_maximal_value = max(np.abs(np.min(matrix)), np.abs(np.max(matrix)))
    im, _ = heatmap(matrix, dc_list, f_list_label, ax=ax, vmin=-bar_maximal_value, vmax=bar_maximal_value,
                    cmap="PuOr", cbarlabel="Velocity(DC, frequency)")
    plt.tight_layout()
    ax.set_title('RIMS: Velocity(DC, frequency) heat-map', fontsize=12, fontweight='bold')
    ax.set_ylabel('Duty Cycle')
    ax.set_xlabel('Ratchet Frequency')

    '''documenting graph'''
    folder = 'Heat maps'
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_name = plot_uid + ' heatmap'
    plt.savefig(folder + r'/' + file_name + '.jpeg')
    print('Graph saved to '+folder+' as ' + file_name)
    plt.close(plot_uid)
    return

def test_decorator():
    headline_panel()
    print("\t-> RIMS TEST MODE <-\nfor normal simulation run program from rims.py")


'''----------------------------------------------------------------------
                            EXECUTION
----------------------------------------------------------------------'''

if __name__ == '__main__':
    test_decorator()
    test_delta_t([150, 200, 250, 300, 350])
