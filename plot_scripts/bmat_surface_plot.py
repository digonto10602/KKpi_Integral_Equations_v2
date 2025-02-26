#This script us used to plot the surface of 
#of the Bmat matrix 

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input


def Bmat_surface(bmat_file, N):
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)
    
    data = np.genfromtxt(bmat_file, unpack=True)

    i_ind = data[0]
    j_ind = data[1]
    re_bmat = data[10]
    im_bmat = data[11]

    xi = np.linspace(i_ind.min(), i_ind.max(), N)
    yi = np.linspace(j_ind.min(), j_ind.max(), N)
    
    Xi, Yi = np.meshgrid(xi, yi) 

    fig, ax = plt.subplots(2,1, figsize = (12, 16))

    Z1 = scipy.interpolate.griddata((i_ind, j_ind), re_bmat, (Xi, Yi), method='linear')
    Z2 = scipy.interpolate.griddata((i_ind, j_ind), im_bmat, (Xi, Yi), method='linear')

    low_limit = -0.1
    hi_limit = 0.1
    color_level = 500
    #ax[0].set_xlim(i_ind.min(), i_ind.max())
    #ax[0].set_ylim(j_ind.min(), j_ind.max())
    ax[0].invert_yaxis()
    ax[0].set_xlim(min(i_ind), max(i_ind))
    ax[0].set_ylim(max(j_ind), min(j_ind)) 
    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    #ax[0].contour(Xi, Yi, Z1, levels=np.linspace(low_limit, hi_limit, num=color_level), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax[0].contourf(Xi, Yi, Z1, levels=np.linspace(low_limit, hi_limit, num=color_level))
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)

    #ax[1].set_xlim(i_ind.min(), i_ind.max())
    #ax[1].set_ylim(j_ind.min(), j_ind.max()) 
    ax[1].invert_yaxis()
    ax[1].set_xlim(min(i_ind), max(i_ind))
    ax[1].set_ylim(max(j_ind), min(j_ind)) 
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].contour(Xi, Yi, Z2, levels=np.linspace(low_limit, hi_limit, num=color_level), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax[0].contourf(Xi, Yi, Z2, levels=np.linspace(low_limit, hi_limit, num=color_level))
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)

    plt.tight_layout()
    outfilename = str(bmat_file).split(".")
    output = outfilename[0] + ".pdf"
    plt.savefig(output)
    #plt.draw()
    #plt.show()
    plt.close()      

if len(sys.argv) > 1:
    filename = sys.argv[1]
    N = int(sys.argv[2])
    bmat_file = filename 
    Bmat_surface(bmat_file, N)
else:
    print("python didn't get proper arguments")
    print("USAGE:<script> <filename> <N>")
    exit() 