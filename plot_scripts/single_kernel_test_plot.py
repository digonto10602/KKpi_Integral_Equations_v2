#This script is used to plot the surface of OPE
#function for G11 G12 and G21. This script will 
#be called from a c++ code that will run this in 
#parallel to generate the surface plots using 
#parallelism. 

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def surface_plot_1(opefile, singularity_file, N):
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)
    
    #The number of points in each axis is hard coded to 250
    #N = 250
    
    DATA1 = np.genfromtxt(opefile, skip_header=1, unpack=True)
    DATA2 = np.genfromtxt(singularity_file, unpack=True) 

    re_p = DATA1[0]
    im_p = DATA1[1]

    Recutplus = DATA2[0]
    Imcutplus = DATA2[1]
    Recutminus = DATA2[2]
    Imcutminus = DATA2[3]


    xi = np.linspace(re_p.min(), re_p.max(), N)
    yi = np.linspace(im_p.min(), im_p.max(), N)
    
    Xi, Yi = np.meshgrid(xi, yi) 

    fig, ax = plt.subplots(figsize = (12, 8))

    imK = DATA1[3]

    Z = scipy.interpolate.griddata((re_p, im_p), imK, (Xi, Yi), method='linear')

    ax.set_xlim(re_p.min(), re_p.max())
    ax.set_ylim(im_p.min(), im_p.max()) 
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    ax.contour(Xi, Yi, Z, levels=np.linspace(-300, 300, num=2000), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax.contourf(Xi, Yi, Z, levels=np.linspace(-30, 30, num=30))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)

    ax.plot(Recutplus, Imcutplus, marker='o',markerfacecolor="red", markersize=1, color="black",linestyle='none', markeredgewidth=2, markeredgecolor='red')
    ax.plot(Recutminus, Imcutminus, marker='o',markerfacecolor="blue", markersize=1, color="black",linestyle='none', markeredgewidth=2, markeredgecolor='blue')
            
    plt.tight_layout()
    outfilename = str(opefile).split(".")
    output = outfilename[0] + ".png"
    #plt.savefig(output)
    plt.show()         
    #plt.close() 

#opefile = 'ope_file_80.dat'

if len(sys.argv) > 1:
    filename = sys.argv[1]
    kern_sing_file = sys.argv[2]
    N = int(sys.argv[3])
else: 
    print("python didn't get proper arguments")
    print("USAGE:<script> <filename> <singularities> <N>")
    exit() 

opefile = filename 
surface_plot_1(opefile, kern_sing_file, N) 