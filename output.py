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
from defines import *


'''----------------------------------------------------------------------
                            IMPLEMENTATION
----------------------------------------------------------------------'''

def create_trace_file(rims_object, ion_subject):
    with open(rims_object.path_for_output + 'simulation trace.csv', newline='', mode='a') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['X0[um]', 'X' + str(ion_subject.points) + '[um]'])


def write_to_trace_file(rims_object, ion_subject):
    with open(rims_object.path_for_output + 'simulation trace.csv', newline='', mode='a') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow([ion_subject.x0, ion_subject.L * ion_subject.arena_count + ion_subject.loc])


def create_log_file(rims_object, ion_subject):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    f = open(rims_object.path_for_output + "RIMS simulation log.txt", "a")
    f.write("RIMS simulation log\n\n")
    f.write("\ttime created: " + str(rims_object.start_time) + "\n")
    f.write("\ttest duration: " + str(datetime.now() - rims_object.start_time) + "\n")

    f.write("\n\tparticles in the system: " + rims_object.ion + "\n")
    f.write("\tdiffusion coefficient: " + str(ion_subject.diffusion) + "[m^2/cm] /10^-9\n")
    f.write("\ttemperature: " + str(TEMPERATURE) + "[k]\n")

    f.write("\nRatchet potential profile\n")
    if rims_object.potential_profile_list[3] == 2:  # sin
        f.write("\tfunction: double sin wave \n")
        f.write("\tamplitude: " + str(
            rims_object.potential_profile_list[1] + rims_object.potential_profile_list[2]) + "[V]\n")

    else:  # saw
        f.write("\tfunction: saw wave \n")
        f.write("\tamplitude: " + str(rims_object.potential_profile_list[1]) + "[V]\n")
    f.write("\twidth: " + str(ion_subject.L) + "[um]\n")
    f.write("\tfrequency: " + str(ion_subject.flash_frequency) + "[Hz]\n")
    f.write("\tperiod: " + str(ion_subject.flash_period) + "[sec]\n")
    f.write("\tduty cycle: " + str(rims_object.dc) + "\n")

    if rims_object.flash_mode == 0:  # ON/OFF
        f.write("\tflash mode: ON/OFF\n")
    else:
        f.write("\tflash mode: + / -\n")

    f.write("\nSimulation settings\n")
    f.write("\tparticles simulated: " + str(rims_object.number_of_simulations) + "\n")
    f.write("\tmeasurements per particle: " + str(ion_subject.points) + "\n")
    f.write("\tintervals (delta_t): " + str(ion_subject.interval) + "e-6[sec]\n")
    f.write("\tfriction coefficient (gamma): " + str(ion_subject.gamma) + "\n")
    f.write("\tresolution: " + str(RESOLUTION) + "\n")
    f.write("\tcurrent: " + str(rims_object.current) + "[A]\n")
    f.close()


def print_log_file(rims_object):
    f = open(rims_object.path_for_output+"RIMS simulation log.txt", "r")
    log = f.read()
    print(log)
    f.close()


def plot_potential_profile(rims_object, x,V,E):
    """
    plots and saves E,V profiles as function of x over 1 cycle
    :param rims_object: simulation instance
    :param x: ratchet cycle x axis
    :param V: potential calculations at x
    :param E: -grad(V) at x
    """
    plot_id = create_unique_id()
    plt.figure(plot_id)
    plt.plot(x, V, label="V(x) potential profile", color=YELLOW)

    plt.suptitle('RIMS: Ratchet potential profile', fontsize=14, fontweight='bold')
    plt.xlabel(r"X [$\mu $m]")
    plt.ylabel(r"V [v]")
    plt.legend(loc="upper left")
    save_plots(rims_object, 'Ratchet potential profile V', plot_id)

    plot_id = create_unique_id()
    plt.figure(plot_id)
    plt.plot(x, E, label=r"E(x) electric field = -$\nabla $V", color=PURPLE)

    plt.suptitle('RIMS: Ratchet potential profile', fontsize=14, fontweight='bold')
    plt.xlabel(r"X [$\mu $m]")
    plt.ylabel(r"E [v/$\mu $m]")
    plt.legend(loc="upper left")
    save_plots(rims_object, 'Ratchet potential profile E', plot_id)

    fig, ax1 = plt.subplots()
    plt.suptitle('RIMS: Ratchet potential profile', fontsize=14, fontweight='bold')

    color = YELLOW
    ax1.set_xlabel(r"X [$\mu $m]")
    ax1.set_ylabel(r"V [v]", color=color)
    ax1.plot(x, V, color=color, label="V(x) potential profile")
    ax1.tick_params(axis='y', labelcolor=color)
    plt.legend(loc="lower right")
    ax2 = ax1.twinx()

    color = PURPLE
    ax2.set_ylabel(r"E [v/$\mu $m]", color=color)
    ax2.plot(x, E, color=color, label=r"E(x) electric field = -$\nabla $V")
    ax2.tick_params(axis='y', labelcolor=color)
    plt.legend(loc="lower left")
    fig.tight_layout()
    plt.savefig(rims_object.path_for_output + 'Ratchet potential profile.jpeg')
    return


def save_plots(rims_object, name, fig_number):
    if not os.path.exists(rims_object.path_for_output):
        os.makedirs(rims_object.path_for_output)
    plt.figure(fig_number)
    plt.savefig(rims_object.path_for_output + name + '.jpeg')
    plt.close(fig_number)
    print(name + ' saved to output plots.')
    return


def create_unique_id():
    now = datetime.now()
    c_time = now.ctime()
    uid = c_time.replace(':','').replace(' ','_')
    return uid
