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
from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list('mycmap', ['teal', 'white'])

def eta_dependent_epsilon(eta, En, qb, sigb, kmax, mi, N):
    pi = np.pi
    termA = (eta*kmax)/(2.0*pi*N)
    termB = (4.0*qb*En*En)/(En*En + mi*mi - sigb)
    return termA*termB 
    
def eta_dependent_epsilon_data(filename, eta, minN, maxN):
    eta_input_file = "eta_epsilon_input_for_" + filename 
    (En, qb, sigb, kmax, mi) = np.genfromtxt(eta_input_file, skip_header=1, unpack=True)
    xval = []
    yval = []

    Narr = np.linspace(minN, maxN, 2000)

    for i in range(0,len(Narr),1):
        Nval1 = Narr[i]
        xval.append(Nval1)
        eps = eta_dependent_epsilon(eta, En, qb, sigb, kmax, mi, Nval1)
        yval.append(eps)
        print(i, eps)
    
    Nval = np.array(xval)
    Epsval = np.array(yval) 

    return Nval, Epsval 

def surface_density_rhophib_N_vs_eps(filename, N):
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    (Nval, eps, delta_rhophib) = np.genfromtxt(filename, unpack=True)

    xi = np.linspace(Nval.min(), Nval.max(), N)
    yi = np.linspace(eps.min(), eps.max(), N)

    Xi, Yi = np.meshgrid(xi, yi) 

    fig, ax = plt.subplots( figsize = (12, 10))

    Z = scipy.interpolate.griddata((Nval, eps), delta_rhophib, (Xi, Yi), method='linear')

    ax.set_ylim([4*10**-4,10**-1])
    ax.set_yscale('log')
    ax.set_xlabel("$N$", fontsize=20)
    ax.set_ylabel("$\epsilon/m_1^2$", fontsize=20)
    #ax.contour(Xi, Yi, Z,levels=np.linspace(0, 1000, num=1000), linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax.contourf(Xi, Yi, Z, levels=np.linspace(0, 5*10**2, num=1000), cmap=cmap, norm='log')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.yaxis.set_ticks_position("both")
    cax.xaxis.set_label_position("bottom")
    cbar.ax.set_yticks([10**-1, 10**0, 10**1, 10**2])
    
    #cax.set_yscale('log')
    cax.tick_params(labelsize=18)

    eta_5_N, eta_5_eps = eta_dependent_epsilon_data(filename, 5, Nval.min(), Nval.max())
    ax.plot(eta_5_N, eta_5_eps, linestyle='solid', linewidth=2, color='goldenrod', label="$\eta=5$")

    eta_15_N, eta_15_eps = eta_dependent_epsilon_data(filename, 15, Nval.min(), Nval.max())
    ax.plot(eta_15_N, eta_15_eps, linestyle='solid', linewidth=2, color='rebeccapurple', label="$\eta=15$")

    eta_25_N, eta_25_eps = eta_dependent_epsilon_data(filename, 25, Nval.min(), Nval.max())
    ax.plot(eta_25_N, eta_25_eps, linestyle='solid', linewidth=2, color='darkred', label="$\eta=25$")

    ax.legend() 
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