''' This script is used to plot
the delta_rho_phib density plot, 
the x axis is the matrix size N, and 
the y axis is the epsilon, we can 
also draw epsilon line based on the eta
value '''

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def surface_density_rhophib_N_vs_eps(filename, N):
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    (Nval, eps, delta_rhophib) = np.genfromtxt(filename, unpack=True)

    xi = np.linspace(Nval.min(), Nval.max(), N)
    yi = np.linspace(eps.min(), eps.max(), N)

    Xi, Yi = np.meshgrid(xi, yi) 

    fig, ax = plt.subplots( figsize = (12, 12))

    Z = scipy.interpolate.griddata((Nval, eps), delta_rhophib, (Xi, Yi), method='linear')

    ax.set_xlabel("$N$", fontsize=20)
    ax.set_ylabel("$\epsilon/m_1^2$", fontsize=20)
    ax.contour(Xi, Yi, Z,levels=np.linspace(0, 100, num=1000), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax.contourf(Xi, Yi, Z, levels=np.linspace(0, 100, num=50))
    divider = make_axes_locatable(ax)
    cax = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)
    plt.tight_layout()

    output = filename + ".png"
    plt.savefig(output)
    plt.close() 

if len(sys.argv) > 1:
    filename = sys.argv[1]
    N = int(sys.argv[2])
else:
    print("usage: script <filename> <row/col size>")
    exit()

surface_density_rhophib_N_vs_eps(filename, N)