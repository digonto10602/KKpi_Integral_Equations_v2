#This function generates the surface plots 
#for M3mat 

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def M3_surface_plots():
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    N = 100

    fig, ax = plt.subplots(2, 2, figsize=(24, 24))
    filename1 = "M3_11_mat.dat"
    filename2 = "M3_12_mat.dat"
    filename3 = "M3_21_mat.dat"
    filename4 = "M3_22_mat.dat"

    (sigp1, sigk1, reM1, imM1) = np.genfromtxt(filename1, unpack=True)

    xi = np.linspace(sigp1.min(), sigp1.max(), N)
    yi = np.linspace(sigk1.min(), sigk1.max(), N)

    Xi, Yi = np.meshgrid(xi, yi) 

    Z = scipy.interpolate.griddata((sigp1, sigk1), reM1, (Xi,Yi), method='linear')

    ax[0,0].tick_params(axis='both', which='major', labelsize=20)
    ax[0,0].tick_params(axis='both', which='minor', labelsize=20)
    #ax[0,0].contour(Xi, Yi, Z, levels=np.linspace(-20000000, 20000000, num=1000), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax[0,0].contourf(Xi, Yi, Z, levels=np.linspace(-3000000, 3000000, num=200))
    divider = make_axes_locatable(ax[0,0])
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)
    
    plt.tight_layout()
    #outfilename = str(opefile).split(".")
    #output = outfilename[0] + ".png"
    #plt.savefig(output)
    plt.show()         
    #plt.close() 

M3_surface_plots() 