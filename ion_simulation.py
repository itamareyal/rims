'''
ion_simulation.py
'''

'''----------------------------------------------------------------------
                                IMPORTS
----------------------------------------------------------------------'''
import numpy as np
import matplotlib.pyplot as plt

from scipy.misc import derivative
import progressbar
from datetime import datetime
import os
import random

from ion_parser import *
#from rims import *





def get_potential_profile(a1, a2, L, x):
    # format of elecrtron experiment, units of Vq
    return a1*np.sin(2*np.pi*x / L) + a2*np.sin(4*np.pi*x / L)

def create_arena_E():
    x=np.linspace(0,0.8) # micrometer
    return plt.plot(x,get_potential_profile(-0.25, -0.05, 0.8, x))


    



class ion:
    def __init__(self, element):
        

        #self.loc = random.uniform(0,np.pi)
        #self.speed = random.uniform(0,np.pi)

        self.mass = 6


        self.intervals = 1
        self.intervals_count = 0
        self.points = 100
        self.inarena = 0
        self.loc_hist = 0

    def create_arena(self,a,b,A):
        x1=np.linspace(0,a) # micrometer
        x2=np.linspace(a,b)
        x=np.linspace(0, a+b)
        f1=A * np.divide(x,a)
        f2=A * np.divide((x-(a+b)),(-b))
        step = 0.5*(np.sign(x-a) +1)
        f=  f1 -step*f1 + step* f2
        plt.plot(x, f)
        plt.xlabel("X [um]")
        plt.ylabel("V [v]")



    #def get_max_x0(self, ):
    


