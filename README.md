# RIMS [Ratchet based Ion Movement Simulator]
## Overview
#### __Second prize winner of Tel-Aviv university EE projects contest of 2021__
#### This project refers to the usage of electric ratchets as a mean of pumping ions selectively. Ratchets are essentially potential profiles that change in space and time at a specific frequency. By exploiting the electric charge of ions, said ratchet can induce ion movement inside a compound. Relaying on different diffusion coefficients of different ions, the ion movement created by the ratchet, can be of opposite direction for different ions. Thus, selectively pumping out only one type of ion.
#### The potential profiles defining the ratchet can vary, changing in amplitude, shape and their time period. In order to predict the effect of different profiles, we created a simulation environment that tests the effect of ratchets on ions. The goal of our project was to build a simulator that calculates the movement of multiple ions under any potential profile inserted to the system. Calculation is ‘Monte Carlo’ based, meaning the initial location of the simulated ion is random.   
#### The product is RIMS (Ratchet based Ions Movement Simulator), a console-based application capable of receiving any potential profile, any physical diffusion coefficient and plotting various graphs of ion movement.
#### Patent number: US20200402782A1
## Usage

#### Executing the main script 'rims.exe' launches the console interface. Creating the simulation is done in 2 steps:

------------------------------------------------------------------------------------------------------
                                Step 1- Configure the system
------------------------------------------------------------------------------------------------------

#### Select the type ions to be simulated by entering their numbers (separated by comma) and hitting ENTER. 

#### Ions can either be selected from the list or manually inserted by their diffusion coefficient.

------------------------------------------------------------------------------------------------------
                                Step 2- Configure the ratchet
------------------------------------------------------------------------------------------------------
#### Select the ratchet potential profile from a pre-prepared csv file saved in 'potential profile sources directory'.
#### In addition, an in-console triangle or sin wave can be defined and used.

#### csv ratchet file must be of the following format: 


X | t0 | t1 | ... | tk 
------------ | ------------- | ------------- | ------------- | -------------
v00 | v01 | v02 | ... | v0n
v10 | v11 | v12 | ... | v1n
v...0 | v...1 | v...2 | ... | v...n
vk0 | vk1 | vk2 | ... | vkn

* x   = length of the potential profile [um]
* ti  = absolute time when vi switches to vi+1 [usec]
* tk  = time period of the whole ratchet marked T
* k+1 = number of different potential profile
* n+1 = number of dots describing each potential profile
* vij = potential at point j of profile i [volt]

#### **NOTE: all potential profiles must be of the same length and same number of points. meaning that every row has exactly n+1 points. n,k > 0**

## Settings
#### Some simulation parameters can be edited in settings.csv file. open in notepad, edit, save and re-launch rims to update the settings
#### More ions & their diffusion coefficients can be permanently added to the step1 via diffusion_coefficients.csv

## Compiled version
#### [Link for pre-compiled versions](https://drive.google.com/drive/folders/1z9EpGJJgFvfM1rzAkMjoKMdQWmRHlILU?usp=sharing)

## Structure
* _rims.py_		    	Top module, hosting class ‘rims’
* _interface.py_		Extraction of data from user and external files
* _ion_simulation.py_	Host of ‘ion’ class, running ion iteration loop
* _current_calc.py_		Calculation of particle speed and current
* _outputs.py_		    Plotting, tracing and logging all outputs
* _Defines.py_		    Hard coded data of physical and simulation constants
* _Setup.py_	    	Enclosing the project into an executable file


## Dependencies
* _numpy_
* _matplotlib_
* _opencv_
* _os_
* _csv_


## Credits
#### Dr. Gideon Segev
#### Eran Weil
#### Itamar Eyal
###### Energy devices lab, Tel-Aviv university
