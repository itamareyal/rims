# RIMS [Ratchet based Ion Movement Simulator]
## Overview
#### __Second prize winner of Tel-Aviv university EE projects contest of 2021__
#### This project refers to the usage of electric ratchets as a mean of pumping ions selectively. Ratchets are essentially potential profiles that change in space and time at a specific frequency. By exploiting the electric charge of ions, said ratchet can induce ion movement inside a compound. Relaying on different diffusion coefficients of different ions, the ion movement created by the ratchet, can be of opposite direction for different ions. Thus, selectively pumping out only one type of ion.
#### The potential profiles defining the ratchet can vary, changing in amplitude, shape and their time period. In order to predict the effect of different profiles, we created a simulation environment that tests the effect of ratchets on ions. The goal of our project was to build a simulator that calculates the movement of multiple ions under any potential profile inserted to the system. Calculation is ‘Monte Carlo’ based, meaning the initial location of the simulated ion is random.   
#### The product is RIMS (Ratchet based Ions Movement Simulator), a console-based application capable of receiving any potential profile, any physical diffusion coefficient and plotting various graphs of ion movement.
#### Patent number: US20200402782A1
## Usage- GUI
#### To run the gui, launch application via gui.py script. Make sure you install all dependencies (listed below under dependencies)
    python gui.py
#### The GUI allows the user to choose a potential profile, ions and run simulations. It also allows to change the simulation settings and view the outputs (including videos).
### RIMS GUI simulation:
<img src="https://github.com/itamareyal/rims/blob/main/system_files/simulation.png?raw=true" width="500">

### RIMS GUI output menu:
<img src="https://github.com/itamareyal/rims/blob/main/system_files/output.png?raw=true" width="500">


## Usage- Terminal (no dependencied needed)

#### Executing the file **rims.exe** launches the console interface. If all dependecies are installed (view dependencies below), you can also run rims.py from the cmd prompt / terminal:
    python rims.py

#### Creating the simulation is done in 2 steps:
`Step 1- Configure the system`

#### Select the type ions to be simulated by entering their numbers (separated by comma) and hitting ENTER. 

#### Ions can either be selected from the list or manually inserted by their diffusion coefficient.


`Step 2- Configure the ratchet`

#### Select the ratchet potential profile from a pre-prepared csv file saved in 'potential profile sources directory'.
#### In addition, an in-console triangle or sin wave can be defined and used.
## Potential profile (Ratchet) Formatting
#### csv ratchet file must be of the following format: 


X | t0 | t1 | ... | tk 
------------ | ------------- | ------------- | ------------- | -------------
v00 | v01 | v02 | ... | v0n
v10 | v11 | v12 | ... | v1n
v...0 | v...1 | v...2 | ... | v...n
vk0 | vk1 | vk2 | ... | vkn

* x   = length of the potential profile [um]
* ti  = absolute time when vi switches to vi+1 [usec]
* tk  = time period of the whole ratchet marked T [usec]
* k+1 = number of different potential profile
* n+1 = number of dots describing each potential profile
* vij = potential at point j of profile i [volt]

#### for example, The below cvs file translates to the following profile:
1|0.6|1| | | | 
------------ | ------------- | ------------- | ------------ | ------------- | ------------- | 
0|0.2|0.35|0.45|0.5|0.3|0

<img src="https://github.com/itamareyal/rims/blob/main/system_files/r1.jpeg?raw=true" width="360">


#### notice how when k=1, RIMS copies the only profile and creates another that is multiplied by -1

#### **NOTE: all potential profiles must be of the same length and same number of points. meaning that every row has exactly n+1 points. dx=x/(n+1) constant for all profiles. n,k > 0**

## Settings
#### Some simulation parameters can be edited in settings.csv file. open in notepad, edit, save and re-launch rims to update the settings
#### More ions & their diffusion coefficients can be permanently added to the step 1 via diffusion_coefficients.csv or the GUI.
#### All of these parameters can also be edited in the GUI.

## Compiled version
#### [Link for pre-compiled versions](https://drive.google.com/drive/folders/1z9EpGJJgFvfM1rzAkMjoKMdQWmRHlILU?usp=sharing)

## Structure
* _rims.py_		    	Top module, hosting class ‘rims’
* _interface.py_		Extraction of data from user and external files
* _ion_simulation.py_	Host of ‘ion’ class, running ion iteration loop
* _current_calc.py_		Calculation of particle speed and current
* _outputs.py_		    Plotting, tracing and logging all outputs
* _defines.py_		    Hard coded data of physical and simulation constants
* _test.py_	    	    Test scripts for multiple simulation runs under different parameters
* _setup.py_	    	Enclosing the project into an executable file
* _gui.py_              Execution with this script as 'main' will launch the RIMS GUI


## Dependencies
* _numpy_
* _matplotlib_
* _opencv_
* _os_
* _csv_
* _PySimpleGui_ (For GUI usage only)

## Outputs
#### All outputs are saved to folder "simulation outputs / [time stamp] [Ion] /"
#### All outputs are also avilable to view via the GUI
* _Distribution histogram of all ions simulated_
  
<img src="https://github.com/itamareyal/rims/blob/main/system_files/Distribution%20histogram.png?raw=true" width="360">

* _Distribution histogram periodic_
  
<img src="https://github.com/itamareyal/rims/blob/main/system_files/Distribution%20histogram%20periodic.jpeg?raw=true" width="360">

* _Distribution histogram infinite_
  
<img src="https://github.com/itamareyal/rims/blob/main/system_files/Distribution%20histogram%20infinite.jpeg?raw=true" width="360">

* _Ratchet potential profiles_
  
<img src="https://github.com/itamareyal/rims/blob/main/system_files/Ratchet%20potential%20profiles.jpeg?raw=true" width="360">

* _Average speed of ions over ratchet cycle_
  
<img src="https://github.com/itamareyal/rims/blob/main/system_files/Average%20speed%20of%20ions%20over%20ratchet%20cycles.jpeg?raw=true" width="360">

* _RIMS simulation summary_
```
RIMS simulation summary

	time created: 2021-10-03 11:22:17.954820
	test duration: 0:00:30.616135
	simulated time: 1.28e-05[sec]

	particles in the system: Lead Pb+2
	diffusion coefficient: 9.45e-06[cm^2/sec]
	temperature: 293[k]

Ratchet potential profile
	width: 1.0[um]
	frequency: 1250000.0[Hz]
	period: 8e-07[sec]

Simulation settings
	number of particles simulated: 2500
	measurements per particle: 240
	time measurement intervals (delta_t): 5.291005291005292e-08[sec]
	friction coefficient (gamma): 2671.829254778837[eVsec/cm^2]
	resolution: 1000 (no. of dx along a single ratchet)
	velocity: -203.63324149682964[cm/sec]
```
* _RIMS simulation log_
```
RIMS simulation log

11:22:17 -	gamma factor calculated 2671.829254778837[eVsec/cm^2]
11:22:17 -	delta t calculated and is 5.291005291005292e-08[sec]
11:22:17 -	electric field generated
11:22:20 -	potential profile plot saved
11:22:20 -	RIMS simulation initialized
11:22:25 -	cycle 0 completed. velocity is -126.53171866273395[cm/sec]
11:22:25 -	cycle count still lower than MIN_MEASUREMENTS_FOR_SS=10
...
11:22:40 -	cycle 9 completed. velocity is -200.33142579255215[cm/sec]
11:22:40 -	steady state reached after 9 cycles
11:22:41 -	cycle 10 completed. velocity is -204.82592543163648[cm/sec]
11:22:41 -	steady state maintained after 9 cycles
...
11:22:48 -	cycle 15 completed. velocity is -198.7913711481517[cm/sec]
11:22:48 -	steady state maintained after 9 cycles
11:22:48 -	data collected for histogram
11:22:48 -	summary file created
11:22:50 -	Distribution histogram infinite saved
11:22:53 -	Distribution histogram periodic saved
11:22:53 -	average speed plot saved
11:22:53 -	simulation finished
11:22:53 -	summary file printed
11:24:11 -	Data added to multiple ions histogram
11:24:17 -	Multiple ions histogram saved
```

* _Simulation trace_
#### Video outputs are saved in "simulation outputs/00_Video outputs" folder
#### Histograms of all ion types on the same plot are saved in "simulation outputs/00_Multiple ions histograms" folder

## Credits
#### Eran Weil
#### Itamar Eyal
#### Dr. Gideon Segev
###### Energy devices lab, Tel-Aviv university

<img src="https://github.com/itamareyal/rims/blob/main/system_files/IMG-3441.jpg?raw=true" width="360">

