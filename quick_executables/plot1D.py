#executable to plot 1-D functions 

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def plot1D_func(filename, ind1, ind2):

    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    DATA = np.genfromtxt(filename, unpack=True)

    x_data = DATA[ind1]
    y_data = DATA[ind2]

    #print(DATA)

    #print("x_data = ",x_data)
    #print("y_data = ",y_data) 

    fig, ax = plt.subplots( figsize = (12, 5) )

    ax.plot(x_data, y_data, marker='o', markersize=20, markeredgewidth=2 )
    plt.tight_layout()
    plt.show()



if len(sys.argv) > 1:
    filename = sys.argv[1]
    ind1 = sys.argv[2]
    ind2 = sys.argv[3]
else:
    print("[USAGE]: python plot1D.py <filename> <x-axis ind> <y-axis ind>")
    exit()

plot1D_func(filename, int(ind1), int(ind2))

