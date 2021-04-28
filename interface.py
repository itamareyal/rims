
from defines import *
import csv
import os

def headline_panel():
    print("-------------------------------------------------------")
    print("     RIMS - Ratchet based Ion Movement Simulator")
    print("-------------------------------------------------------")

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

    print("\nEnter ratchet function:\n\t0)Load from csv sheet\n\t1)Saw wave\n\t2)Double Sin\n")
    ratchet_number = input_check_int("Ratchet function number = ", [0, 1, 2])
    if ratchet_number != 0:
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
            x = np.linspace(0, a + b, num=RESOLUTION)
            f1 = A * np.divide(x, a)
            f2 = A * np.divide((x - (a + b)), (-b))
            step = np.heaviside(x - a, 1)
            pos = f1 - step * f1 + step * f2

        elif ratchet_number == 2:
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

            x = np.linspace(0, L, num=RESOLUTION)
            pos = a1 * np.sin(2 * np.pi * x / L) + a2 * np.sin(4 * np.pi * x / L)
        neg = np.multiply(pos, -ALPHA)
        potential_mat = np.vstack((pos, neg))

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
        T = 1 / flash_frequency
        t_vec = np.array([dc * T, T])

    else:                                                   # data from csv
        L, t_vec, potential_mat = select_csv_file()
        x = np.linspace(0, L, num=potential_mat.shape[1])


    potential_profile = [L, x, t_vec, potential_mat]

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
    return ion_selection, potential_profile, output_selection

def ion_selection_panel():
    print("-------------------------------------------------------")
    print("             Step 1- Configure the system")
    print("-------------------------------------------------------")

    print("\nSelect an ion to be simulated from the following list:")
    print("\t1) Lead Pb+2")
    print("\t2) Potassium K+")
    print("\t3) Calcium Ca+2")
    print("\t4) Sodium Na+")
    print("\t5) Electrons in Silicon")
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


def select_csv_file():
    folder = r'potential profile sources/'
    if not os.path.exists(folder):
        os.makedirs(folder)

    items = os.listdir(folder)
    i = 1
    valid_files = []
    print("\nPotential profiles csv/txt files found in folder '"+folder+"':")
    for entry in items:
        if os.path.isfile(os.path.join(folder, entry)) and entry.endswith('csv') or entry.endswith('txt'):
            print(str(i)+') ' + entry)
            i+=1
            valid_files.append(entry)
    if len(valid_files)==0:
        print("\nNo .csv or .txt files found in "+folder)
        exit()
    else:
        print("Select the number of the profile you wish load:")
        file = valid_files[input_check_int("file number: ",range(1,len(valid_files)+1))-1]
    scalar_x, vec_t, mat_v = load_data_from_csv(folder + file)
    return scalar_x, vec_t, mat_v


def load_data_from_csv(csv_file_path):
    """
    :param csv_file_path - path to a ratchet profile setting file set as follows
    x, t0, t1... tn
    v00, v01, v02... v0n
    v10, v11, v12... v1n
    ...
    vk0, vk1, vk2... vkn

    x   = length of the potential profile [um]
    ti  = absolute time when vi switches to vi+1 [SEC]
    tn  = time period of the whole ratchet marked T
    k+1 = number of different potential profile
    vi  = vector describing the shape of potential profile i [volt]

    NOTE: all potential profiles must be of the same length and same number of points
    """
    vec_header = np.loadtxt(csv_file_path, max_rows=1, delimiter=',')   # holds the first data row
    scalar_x = vec_header[0] * pow(10, -4)                              # width of the profile
    vec_t = vec_header[1:]                                              # timings vector
    mat_v = np.loadtxt(csv_file_path, skiprows=1, delimiter=',')        # potential profiles
    if mat_v.ndim==1:
        print("only 1 profile was detected. adding a second profile such that V2= -ALPHA*V1")
        mat_v = np.vstack((mat_v, -ALPHA * mat_v))

    return scalar_x, vec_t, mat_v
