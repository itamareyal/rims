import PySimpleGUI as sg
from PySimpleGUI.PySimpleGUI import No
import numpy
from rims import *


'''----------------------------------------------------------------------
                        SIZES & FONTS
----------------------------------------------------------------------'''
BUTTON_SIZE = (15, 1.1)
PROGRESS_BAR_SIZE = (10, 0.5)
OUTPUT_SIZE = (100, 6)
LISTBOX_SIZE = (80, 3)
FILELIST_SIZE = (25, 12)
TEXT_OUTPUT_SIZE = (600, None)
MENU_FONT = 'Any 15'
TAB_FONT = 'Any 25'
BUTTON_FONT = 'Any 15'

'''----------------------------------------------------------------------
                            FUNCTIONS
----------------------------------------------------------------------'''

def play_video(video_file):
    cap = cv2.VideoCapture(video_file)

    fps = 1000
    frames = []

    #Play the video once and store the frames in an array
    while cap.isOpened():
        _ret, frame = cap.read()        
        if frame is None:
            break
        frames.append(frame)
    cap.release()
    end = False
    while True:
        if end:
            break
        for frame in frames:
            cv2.imshow("ESC to close",frame)
            k = cv2.waitKey(int(1000)) & 0xFF

            #End with ESC
            if k == 27:
                end = True
                break
    cap.release()
    cv2.destroyAllWindows()

def preview_ratchet(l_s, x_v, t_v, v_m, filename):
    "Creates an image preview to be displayed in gui"
    number_of_profiles = v_m.shape[0]
    x = [dx * pow(10, 4) for dx in x_v]
    dx = l_s/load_one_setting(settings_filename,'RESOLUTION')
    
    fig, axs = plt.subplots(number_of_profiles)
    plt.suptitle('Preview: '+filename, fontsize=10, fontweight='bold')

    for i_profile in range(number_of_profiles):
        y = [dy for dy in v_m[i_profile]]
        axs[i_profile].plot(x, y, color=YELLOW, label="V(x) potential profile")
        axs[i_profile].tick_params(axis='y', labelcolor=YELLOW)

    text_kwargs = dict(fontsize=10, color=YELLOW, fontweight='bold')
    fig.text(0.1, 0.91, 'V(x) potential profile', text_kwargs)
    axs[number_of_profiles-1].set_xlabel(r"X [$\mu $m]")
    axs[int(number_of_profiles/2)].set_ylabel(r"V [v]")

    if number_of_profiles < 4:
        fig.tight_layout()
    if filename.endswith('.csv'):
        filename = filename[:-4]
    plt.savefig(filename + '.jpeg')
    plt.close()

def load_first_preview(folder):
    items = os.listdir(folder)
    for item in items:
        if (item.endswith('csv') or item.endswith('txt')):
            scalar_x, vec_t, mat_v = load_data_from_csv(folder+r'/'+item)
            vec_x = np.linspace(0, scalar_x, num=mat_v.shape[1])
            preview_ratchet(scalar_x, vec_x, vec_t, mat_v, item)
            return folder+r'/'+item[:-4]+'.jpeg'
            
def r_print(line):
    if not __name__ == '__main__':
        print(line)
    else:
        sg.Print(line)

def display_folder_csv(folder):
    try:
        file_list = os.listdir(folder)
    except:
        file_list = []

    fnames = [
        f
        for f in file_list
        if os.path.isfile(os.path.join(folder, f))
        and f.lower().endswith((".csv"))
    ]
    return fnames

def display_folder_all(folder):
    try:
        file_list = os.listdir(folder)
    except:
        file_list = []
    fnames = [
        f
        for f in file_list
        
        if (os.path.isfile(os.path.join(folder, f))
        and (f.lower().endswith((".csv")) 
        or f.lower().endswith((".png")) 
        or f.lower().endswith((".jpeg"))
        or f.lower().endswith((".avi"))
        or f.lower().endswith((".txt")))
        or os.path.isdir(os.path.join(folder, f.lower())))
        and not f.lower().startswith('cycle_')
    ]
    fnames = sorted(fnames, key=sort_output_folders)
    return fnames

potential_profiles_folder = "potential profile sources"
histogram_folder = "simulation outputs"

'''----------------------------------------------------------------------
                            LAYOUTS
----------------------------------------------------------------------'''

# --------------------------- Potential profile tab --------------------------- #
file_list_column = [
    [   sg.In( enable_events=True, key="-FOLDER-", default_text=potential_profiles_folder, size=(20,1), font=MENU_FONT),
        sg.FolderBrowse(initial_folder=potential_profiles_folder,button_text="Browse", font=MENU_FONT),
    ],
    [   sg.Listbox(values=display_folder_csv(potential_profiles_folder), enable_events=True, key="-FILE LIST-",font=MENU_FONT)
    ],
]

image_viewer_column = [
    [   sg.Text("Choose a ratchet from list on left:",font=MENU_FONT)],
    [   sg.Text( key="-TOUT-",font=MENU_FONT)],
    [   sg.Image(key="-IMAGE-",filename='system_files/blank_preview.jpeg')],
]

layout_potential_profile = [
    [
        sg.Column(file_list_column, size=(10,1)),
        sg.VSeperator(),
        sg.Column(image_viewer_column, size=(10,1)),
    ]
]

# --------------------------- Ion selection tab --------------------------- #
layout_ions = [
    [   sg.Text("Choose ion/s from the list and click 'Select':",font=MENU_FONT)],
    [   sg.Listbox(values=diffusion_coefficient_dict.items(),size=LISTBOX_SIZE, select_mode='extended', key="lstbox_ions", font=MENU_FONT)],
    [   sg.Button('Select', size=BUTTON_SIZE, font=BUTTON_FONT),sg.Button('Reset',size=BUTTON_SIZE, font=BUTTON_FONT),
    sg.Button('Add more...',size=BUTTON_SIZE, font=BUTTON_FONT), sg.Button('Remove',size=BUTTON_SIZE, font=BUTTON_FONT)],
    [   sg.Text("Ions selected:",font=MENU_FONT)],
    [   sg.Listbox(values=[], key="lstbox_ions_selected",font=MENU_FONT,size=LISTBOX_SIZE)]
]

# --------------------------- Settings tab --------------------------- #
layout_settings = [
    [   sg.Text("Please don't edit settings while the simulation is running",text_color='darkred',font=MENU_FONT)],
    [   sg.Checkbox("Add video output",tooltip='Generate a video of ion distribution over time.',font=MENU_FONT, key='enable_video',enable_events=True,default=load_one_setting(settings_filename,'ENABLE_VIDEO'))],
    [   sg.Checkbox("Manually set time interval", font=MENU_FONT, key='overwrite_delta_t',enable_events=True,default=load_one_setting(settings_filename,'OVERWRITE_DELTA_T')),
    sg.In( enable_events=True, key="delta_t", font=MENU_FONT, default_text=load_one_setting(settings_filename,'delta_t')),
    sg.Text("[usec]")],
    [   sg.Checkbox("Add full trace file", font=MENU_FONT, tooltip='Create csv file with all particle locations at all time intervals',key='create_trace',enable_events=True,default=load_one_setting(settings_filename,'CREATE_TRACE'))],
    [   sg.Text("Particles simulated",font=MENU_FONT, tooltip='The number of particles of each ion type to be simulated'),
    sg.Slider(range=(1, 20),text_color='white',orientation='h', size=(34, 20), default_value=6, tick_interval=2, enable_events=True, key='particles_simulated', font=MENU_FONT),sg.Text(load_one_setting(settings_filename,'particles_simulated'),key='particles_simulated_text', font=MENU_FONT)],
    [   sg.Text("Steady state percent margin", font=MENU_FONT,tooltip='the error margin allowed for steady state'),sg.In(enable_events=True,default_text=load_one_setting(settings_filename,'steady_state_percent_margin') ,key="steady_state_percent_margin", font=MENU_FONT),sg.Text("%",font=MENU_FONT)],
    [   sg.Text("Maximal number of ratchet cycles", font=MENU_FONT),sg.In(enable_events=True,default_text=load_one_setting(settings_filename,'MAX_CYCLES') ,key="max_cycles", font=MENU_FONT)],
    [   sg.Text("Resolution", font=MENU_FONT, tooltip='Resolution of histogram. 1000 is recommended'),sg.In(enable_events=True,default_text=load_one_setting(settings_filename,'resolution') ,key="resolution", font=MENU_FONT)],
    [   sg.Text("Ratchets in display", font=MENU_FONT, tooltip='The number of ratchets in space displayed in the periodic histogram. has no impact on results.'),sg.In(enable_events=True,default_text=load_one_setting(settings_filename,'ratchets_in_system') ,key="ratchets_in_system", font=MENU_FONT)],
]
# --------------------------- Outputs tab --------------------------- #
output_list_column = [
    [   sg.In( enable_events=True, key="-RESULTS_FOLDER-", default_text=histogram_folder, size=(20,1), font=MENU_FONT),
        sg.FolderBrowse(initial_folder=histogram_folder, button_text="Browse", font=MENU_FONT),
    ],
    [
        sg.Button('Simulation outputs folder', font=MENU_FONT, size=(20,1))
    ],
    [   sg.Listbox(values=display_folder_all(histogram_folder), enable_events=True, key="-RESULTS_LIST-",font=MENU_FONT)
    ],
]

output_tab_img = [[   sg.Image(key="-RESULTS_IMAGE-",filename='system_files/blank_preview.jpeg')]]
output_tab_txt = [[   sg.Multiline(key="-RESULTS_MLINE-", font=MENU_FONT,size=TEXT_OUTPUT_SIZE)]]

output_tabgroup = [[sg.TabGroup([
    [sg.Tab('img', output_tab_img, key='tab_img'),
    sg.Tab('txt', output_tab_txt, key='tab_txt'),]
])]]

layout_output = [
    [
        sg.Column(output_list_column, size=(10,1)),
        sg.VSeperator(),
        sg.Column(output_tabgroup, size=(10,1)),
    ]
]

# --------------------------- Full layout --------------------------- #
tabgroup = sg.TabGroup([
    [   sg.Tab('Potential profile', layout_potential_profile,font=TAB_FONT),
        sg.Tab('Ions selection',layout_ions,font=TAB_FONT),
        sg.Tab('Settings', layout_settings,font=TAB_FONT),
        sg.Tab('Outputs', layout_output,font=TAB_FONT)
    ]],font=TAB_FONT,tab_location='lefttop',title_color='white', background_color='#8B8B83', selected_title_color='#607B8B',
    ) 
        
layout = [
    [   tabgroup],
    [   sg.Button('Run Simulation',size=BUTTON_SIZE, font=BUTTON_FONT), sg.Text("", key='Status', font='Any 20', text_color='darkred')],
    [   sg.ProgressBar(1000, orientation='h', size=PROGRESS_BAR_SIZE, key='progressbar',visible=False)],
    [   sg.Output(key='output', visible=False, size=OUTPUT_SIZE)]
    ]

window = sg.Window("RIMS - Ratchet based Ion Movement Simulator", layout)
selected_potential_profile = None
ratchetfile = None
selected_ions = {}

# --------------------------- Main GUI loop --------------------------- #
while True:
    event, values = window.read()
    if event == "Exit" or event == sg.WIN_CLOSED:
        break

    # --------------------------- Potential profile events --------------------------- #
    if event == "-FOLDER-":
        folder = values["-FOLDER-"]
        window["-FILE LIST-"].update(display_folder_csv(folder))
    elif event == "-FILE LIST-":
        try:
            ratchetfile = values["-FILE LIST-"][0]
            filename = os.path.join(
                values["-FOLDER-"], values["-FILE LIST-"][0]
            )
            scalar_x, vec_t, mat_v = load_data_from_csv(filename)
            vec_x = np.linspace(0, scalar_x, num=mat_v.shape[1])
            preview_ratchet(scalar_x, vec_x, vec_t, mat_v, filename)
            window["-TOUT-"].update("Presenting: "+filename)
            window["-IMAGE-"].update(filename=filename[:-4]+'.jpeg')
            selected_potential_profile = [scalar_x, vec_x, vec_t, mat_v]
        except:
            pass

    # --------------------------- Ion selection events --------------------------- #
    elif event == 'Select':
        for ion_box in values["lstbox_ions"]:
            selected_ions[ion_box[0]] = ion_box[1]
        window["lstbox_ions_selected"].update(values=selected_ions.items())
    elif event == 'Reset':
        selected_ions = {}
        window["lstbox_ions_selected"].update(values=selected_ions.items())
    elif event == 'Add more...':
        new_ion_tup = sg.popup_get_text('Insert new ion in the format: ion type, diffusion coefficient [cm^2/sec]',title='Insert new ion',default_text='ion type, diffusion coefficient',size=(50,1),font=MENU_FONT)
        if new_ion_tup in ['ion type, diffusion coefficient',None]:
            continue
        if add_ion_to_dict(ion_file, new_ion_tup):
            diffusion_coefficient_dict = load_settings(ion_file)
            window["lstbox_ions"].update(values=diffusion_coefficient_dict.items())
        else:
            sg.popup('Invalid ion input. Enter ion type and D seperated by comma',title='Notification',font=MENU_FONT)
    elif event=='Remove':
        if len(values["lstbox_ions"])==0:
            sg.popup('Select an entry to remove then click Remove',font=MENU_FONT,title='Notification')
        elif len(values["lstbox_ions"]) >= len(diffusion_coefficient_dict.items()):
            sg.popup('Cannot delete ALL the ions on the list',font=MENU_FONT,title='Notification')
        else:
            for ion_box in values["lstbox_ions"]:
                remove_ion_from_dict(ion_file, ion_box)
                diffusion_coefficient_dict = load_settings(ion_file)
                window["lstbox_ions"].update(values=diffusion_coefficient_dict.items())

    # --------------------------- Settings events --------------------------- #
    elif event == 'enable_video':
        settings = edit_settings(settings_filename,'enable_video',values['enable_video'])
    elif event == 'particles_simulated':
        window['particles_simulated_text'].update(str(values['particles_simulated']*500))
        settings = edit_settings(settings_filename,'particles_simulated',str(values['particles_simulated']*500))
    elif event == 'steady_state_percent_margin':
        settings = edit_settings(settings_filename,'steady_state_percent_margin',str(values['steady_state_percent_margin']))
    elif event == 'max_cycles':
        settings = edit_settings(settings_filename,'max_cycles',str(values['max_cycles']))
    elif event == 'ratchets_in_system':
        settings = edit_settings(settings_filename,'ratchets_in_system',str(values['ratchets_in_system']))
    elif event == 'overwrite_delta_t':
        settings = edit_settings(settings_filename,'overwrite_delta_t',str(values['overwrite_delta_t']))
    elif event == 'delta_t':
        settings = edit_settings(settings_filename,'delta_t',str(values['delta_t']))
    
    # --------------------------- Simulation events --------------------------- #
    elif event == 'Run Simulation':
        if selected_potential_profile != None and selected_ions!={}:
            window['output'].update(visible=True)
            window['progressbar'].update(visible=True)
            window['Status'].update("Simulation in progress...")
            window.Refresh()
            window.Refresh()
            execution(selected_ions, selected_potential_profile, window)
            window['Status'].update("Simulation finished!")
        else:
            window['Status'].update("Select both Potential profile and ion/s")

    # --------------------------- outputs events --------------------------- #
    elif event == "-RESULTS_FOLDER-":
        folder = values["-RESULTS_FOLDER-"]
        window["-RESULTS_LIST-"].update(display_folder_all(folder))
    elif event == "-RESULTS_LIST-":
        try:
            output_file = values["-RESULTS_LIST-"][0]
            filename = os.path.join(
                values["-RESULTS_FOLDER-"], values["-RESULTS_LIST-"][0]
            )
            if filename.lower().endswith('.txt') or filename.lower().endswith('.csv'):
                with open(filename) as txt_file:
                    lines = txt_file.read()

                window["-RESULTS_MLINE-"].update(lines)
                window['tab_txt'].select()
                #window["-RESULTS_MLINE-"].update(visible=True)
                #txt_file.close()
                #window["-RESULTS_IMAGE-"].update(disabled=True)
            elif filename.lower().endswith('.jpeg') or filename.lower().endswith('.png'):
                window["-RESULTS_IMAGE-"].update(filename=filename)
                window['tab_img'].select()

            elif filename.lower().endswith('.avi'):
                play_video(filename)

            elif os.path.isdir(filename.lower()):
                folder = filename
                window["-RESULTS_FOLDER-"].update(folder)
                window["-RESULTS_LIST-"].update(display_folder_all(folder))
        except:
            pass
    elif event == 'Simulation outputs folder':
        window["-RESULTS_FOLDER-"].update(histogram_folder)
        window["-RESULTS_LIST-"].update(display_folder_all(histogram_folder))


window.close()

