"""
output.py

Creates and writes data to log, trace and plot
"""

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''

import matplotlib.pyplot as plt
from datetime import datetime
import os
import csv
import sys
from defines import *


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''


def create_trace_file(rims_object):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    with open(rims_object.path_for_output + 'simulation trace.csv', newline='', mode='a') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['X0[cm]', 'X' + str(POINTS) + '[cm]'])


def create_test_csv(rims_object):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    with open(rims_object.path_for_output + 'test_profiles.csv', newline='', mode='a') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        t = [0.8, 0.001, 0.002, 0.004]
        writer.writerow(t)
        writer.writerow(rims_object.potential_profile)
        writer.writerow(-v for v in rims_object.electric_field)
        writer.writerow(v+1 for v in rims_object.electric_field)
    csv_file.close()
    print("test csv printed")


def write_to_trace_file(rims_object, ion_subject):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    with open(rims_object.path_for_output + 'simulation trace.csv', newline='', mode='a') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow([ion_subject.x0, ion_subject.L * ion_subject.arena_count + ion_subject.loc])


def create_log_file(rims_object):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    f = open(rims_object.path_for_output + "RIMS simulation log.txt", "a")
    f.write("RIMS simulation log\n\n")
    f.write("\ttime created: " + str(rims_object.start_time) + "\n")
    f.write("\ttest duration: " + str(datetime.now() - rims_object.start_time) + "\n")
    f.write("\tsimulated time: " + str(rims_object.cycles_count * rims_object.flash_period) + "[sec]\n")

    f.write("\n\tparticles in the system: " + rims_object.ion + "\n")
    f.write("\tdiffusion coefficient: " + str(rims_object.diffusion) + "[cm^2/sec]\n")
    f.write("\ttemperature: " + str(TEMPERATURE) + "[k]\n")

    f.write("\nRatchet potential profile\n")
    f.write("\twidth: " + str(rims_object.L * pow(10, 4)) + "[um]\n")
    f.write("\tfrequency: " + str(rims_object.flash_frequency) + "[Hz]\n")
    f.write("\tperiod: " + str(rims_object.flash_period) + "[sec]\n")

    if len(rims_object.time_vec) == 2:                              # for 2 state ratchet
        f.write("\tduty cycle: " + str(rims_object.time_vec[0] * rims_object.flash_frequency) + "\n")

    f.write("\nSimulation settings\n")
    f.write("\tparticles simulated: " + str(rims_object.number_of_simulations) + "\n")
    f.write("\tmeasurements per particle: " + str(rims_object.cycles_count * rims_object.intervals_in_period) + "\n")
    f.write("\tintervals (delta_t): " + str(rims_object.interval) + "[sec]\n")
    f.write("\tfriction coefficient (gamma): " + str(rims_object.gamma) + "\n")
    f.write("\tresolution: " + str(rims_object.resolution) + " (no. of dx along a single ratchet)\n")
    f.write("\tvelocity: " + str(rims_object.velocity) + "[cm/sec]\n")
    f.close()


def print_log_file(rims_object):
    f = open(rims_object.path_for_output+"RIMS simulation log.txt", "r")
    log = f.read()
    print(log)
    f.close()

def plot_potential_profile(rims_object):
    """
    plots and saves E,V profiles as function of x over 1 cycle
    :param rims_object: simulation instance
    """
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    number_of_profiles = rims_object.potential_profile_mat.shape[0]
    x = [dx * pow(10, 4) for dx in rims_object.x_space_vec]
    fig, axs = plt.subplots(number_of_profiles)
    plt.suptitle('RIMS: Ratchet potential profiles', fontsize=14, fontweight='bold')

    for i_profile in range(number_of_profiles):
        y = [dy for dy in rims_object.potential_profile_mat[i_profile]]
        axs[i_profile].plot(x, y, color=YELLOW, label="V(x) potential profile")
        axs[i_profile].tick_params(axis='y', labelcolor=YELLOW)
        ax2 = axs[i_profile].twinx()
        y = [dy * pow(10, -4) for dy in rims_object.electric_field_mat[i_profile]]
        ax2.plot(x, y, color=PURPLE, label=r"E(x) electric field = -$\nabla $V")
        ax2.tick_params(axis='y', labelcolor=PURPLE)
    text_kwargs = dict(fontsize=10, color=YELLOW, fontweight='bold')
    fig.text(0.1, 0.91, 'V(x) potential profile [v]', text_kwargs)
    text_kwargs = dict(fontsize=10, color=PURPLE, fontweight='bold')
    fig.text(0.5, 0.91, r"E(x) electric field = -$\nabla $V [v/$\mu $m]", text_kwargs)
    axs[number_of_profiles-1].set_xlabel(r"X [$\mu $m]")

    if number_of_profiles < 4:
        fig.tight_layout()
    plt.savefig(rims_object.path_for_output + 'Ratchet potential profiles.jpeg')
    plt.close()
    return

def plot_average_speed_of_ions(rims_object, v_plot_list):
    unique_id = create_unique_id()
    plt.figure(unique_id)
    x_axis = [cycle + 1 for cycle in range(len(v_plot_list))]
    plt.plot(x_axis, v_plot_list)
    plt.xlabel(r"Ratchet Cycle")
    plt.ylabel(r"Particle Velocity [cm/sec]")
    plt.suptitle("RIMS: Average speed of ions over ratchet cycles", fontsize=14, fontweight='bold')
    save_plots(rims_object, 'Average speed of ions over ratchet cycles',unique_id)
    return

def plot_distribution_over_x_histogram(rims_object, x_plot_list):
    """
    :param x_plot_list: array of final absolute locations of ions
    :param rims_object: simulation object
    """
    '''Infinite system'''
    plot_id = create_unique_id()
    plt.figure(plot_id)
    weights = np.ones_like(x_plot_list) / float(len(x_plot_list))
    x_results_um = [x * np.power(10, 4) for x in x_plot_list]
    plt.hist(x_results_um, weights=weights, bins=rims_object.resolution, label=r"X [$\mu $m]")
    plt.ylabel('Density')
    plt.xlabel(r'X [$\mu $m]')
    plt.title(r"RIMS: Histogram of distribution in infinite system: $\rho $(x)", fontsize=14, fontweight='bold')
    save_plots(rims_object, 'Distribution x axis histogram infinite', plot_id)
    plt.close(plot_id)

    '''keeping final location in range [0, num of ratchets times arena size] to view steady state'''
    plot_id = create_unique_id()
    plt.figure(plot_id)
    x_periodic_system = []
    for x in x_plot_list:
        x_periodic = x
        while x_periodic > RATCHETS_IN_SYSTEM * rims_object.L:
            x_periodic -= RATCHETS_IN_SYSTEM * rims_object.L
        while x_periodic < 0:
            x_periodic += RATCHETS_IN_SYSTEM * rims_object.L
        x_periodic_system.append(x_periodic)
    x_results_um = [x * np.power(10, 4) for x in x_periodic_system]
    plt.hist(x_results_um, weights=weights, bins=rims_object.resolution, label=r"X [$\mu $m]")
    plt.ylabel('Density')
    plt.xlabel(r'X [$\mu $m]')
    plt.title(r"RIMS: Histogram of distribution in periodic system: $\rho $(x)", fontsize=14, fontweight='bold')
    save_plots(rims_object, 'Distribution x axis histogram periodic', plot_id)
    return


def save_plots(rims_object, name, fig_number):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    plt.figure(fig_number)
    plt.savefig(rims_object.path_for_output + name + '.jpeg')
    plt.close(fig_number)
    print('\n' + name + ' saved to output plots.')
    return


def create_unique_id():
    now = datetime.now()
    c_time = now.ctime()
    uid = c_time.replace(':', '').replace(' ', '_')
    return uid


def percentage_progress(n, N):
    progress = int(n * 100) / int(N)
    sys.stdout.write("\r%d%%" % progress)
    sys.stdout.flush()


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar