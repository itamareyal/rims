
from defines import *

def headline_panel():
    print("-------------------------------------------------------")
    print("     RIMS - Ratchet based Ion Movement Simulator")
    print("-------------------------------------------------------")

    input("\npress ENTER to begin...\n")
    return

def execution_rerun_panel():
    print("\n-------------------------------------------------------")
    print("                 Simulation over")
    print("-------------------------------------------------------")
    rerun = input("\nPress y and ENTER to run an new simulation, otherwise press ENTER...\n")
    if rerun == 'y':
        return True
    return False

def extract_data_from_interface(number_selection):
    ion_selection = ION_LIST[number_selection - 1]

    print("\nIon selected: " + ion_selection)

    print("-------------------------------------------------------")
    print("             Step 2- Configure the ratchet")
    print("-------------------------------------------------------")

    print("\nEnter ratchet function:\n\t1)Saw wave\n\t2)Double Sin\n")
    ratchet_number = input_check_int("Ratchet function number = ", [1, 2])
    if ratchet_number == 1:

        print("Please describe the potential profile applied on the system.\n")
        print("                    _          ")
        print("        / |\         |         ")
        print("      /   | \        | A[v]    ")
        print("    /     |  \       |         ")
        print("  /       |   \     _|         ")
        print("                               ")
        print("  \______/\____/               ")
        print("    a[um]   b[um]            \n")

        a = input_check_float("\ta[um] = ") * pow(10, -4)
        b = input_check_float("\tb[um] = ") * pow(10, -4)
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
        L = input_check_float("\tL[um] = ") * pow(10, -4)
        a1 = input_check_float("\ta1[v] = ")
        a2 = input_check_float("\ta2[v] = ")
        potential_profile = [L, a1, a2, ratchet_number]

    print("\nEnter ratchet flashing frequency in Hz: (can use factors of K,M,G)")
    flash_frequency = -1
    tries = 0
    while flash_frequency <= 0:
        if tries > 0:
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
    flash_mode = FLASHING_MODES[flash_number - 1]

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

    print("Initial location of every ion is uniformly randomized over [0," + str(L * pow(10, 4)) + "um]\n")
    return ion_selection, potential_profile, flash_frequency, flash_mode, dc, output_selection

def ion_selection_panel():
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
    number_selection = input_check_int("Enter your selection:", range(1, 7))
    return number_selection


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
