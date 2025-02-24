#This scripts plots functions
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

folder1 = "/home/digonto/Codes/Practical_Lattice_v2/OUTPUTS/N_dependence/"

def N_dependence_plots_1():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = folder1 + "mphib_BFmethod_vs_N_En_eta_25_1_forpaper.dat"
    filename2 = folder1 + "mphib_BFmethod_vs_N_En_eta_25_2_forpaper.dat"
    filename3 = folder1 + "Mphib_SA_method_vs_N_eta_25_1_1_forpaper.dat"
    filename4 = folder1 + "Mphib_SA_method_vs_N_eta_25_1_2_forpaper.dat"

    (N1, diff1, Re_rhoM_1, Im_rhoM_1) = np.genfromtxt(filename1, unpack=True)
    (N2, diff2, Re_rhoM_2, Im_rhoM_2) = np.genfromtxt(filename2, unpack=True)
    (N3, diff3, Re_rhoM_3, Im_rhoM_3) = np.genfromtxt(filename3, unpack=True)
    (N4, diff4, Re_rhoM_4, Im_rhoM_4) = np.genfromtxt(filename4, unpack=True)

    fig, ax = plt.subplots(2, 2, figsize = (24, 12))

    ax[0,0].set_title("BF Method")
    ax[0,0].set_ylim([-0.02,1.0])
    ax[0,0].plot(N1, Re_rhoM_1, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,0].plot(N1, Im_rhoM_1, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,0].plot(N2, Re_rhoM_2, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,0].plot(N2, Im_rhoM_2, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[0,0].set_ylabel("$|\\rho_{\\varphi b} \\mathcal{M}_{\\varphi b}|$")

    ax[1,0].set_ylim([10,40])
    ax[1,0].plot(N1, diff1, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].plot(N2, diff2, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].set_ylabel("$\\Delta \\rho_{\\varphi b}$")
    ax[1,0].set_xlabel("$N$")

    ax[0,1].set_title("SA Method")
    ax[0,1].set_ylim([-0.02,1.0])
    ax[0,1].plot(N3, Re_rhoM_3, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,1].plot(N3, Im_rhoM_3, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,1].plot(N4, Re_rhoM_4, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,1].plot(N4, Im_rhoM_4, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[1,1].set_ylim([10,40])
    ax[1,1].plot(N3, diff3, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,1].plot(N4, diff4, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    #ax[1,1].set_ylabel("$difference$")
    ax[1,1].set_xlabel("$N$")

    #ax[0,0].set_xscale('log')
    #ax[0,0].set_yscale('log')
    #ax[0,1].set_xscale('log')
    #ax[0,1].set_yscale('log')
    #ax[1,0].set_xscale('log')
    #ax[1,0].set_yscale('log')
    #ax[1,1].set_xscale('log')
    #ax[1,1].set_yscale('log')

    ax[0,0].legend()
    ax[0,1].legend()
    #ax[0,0].legend()
    #ax[0,0].legend()
    plt.tight_layout()
    output = "difference_between_bf_and_sa_method.pdf"
    plt.savefig(output)
    plt.show() 


def N_dependence_plots_2():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = folder1 + "mphib_BFmethod_vs_N_En_eta_25_3_forpaper.dat"
    #filename2 = folder1 + "mphib_BFmethod_vs_N_En_eta_25_2_forpaper.dat"
    filename3 = folder1 + "Mphib_SA_method_vs_N_eta_25_1_3_forpaper.dat"
    #filename4 = folder1 + "Mphib_SA_method_vs_N_eta_25_1_2_forpaper.dat"

    (N1, diff1, Re_rhoM_1, Im_rhoM_1) = np.genfromtxt(filename1, unpack=True)
    #(N2, diff2, Re_rhoM_2, Im_rhoM_2) = np.genfromtxt(filename2, unpack=True)
    (N3, diff3, Re_rhoM_3, Im_rhoM_3) = np.genfromtxt(filename3, unpack=True)
    #(N4, diff4, Re_rhoM_4, Im_rhoM_4) = np.genfromtxt(filename4, unpack=True)

    fig, ax = plt.subplots(2, 2, figsize = (24, 12))

    ax[0,0].set_title("BF Method")
    ax[0,0].set_ylim([-0.1,1.0])
    #ax[0,0].plot(N1, Re_rhoM_1, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    #ax[0,0].plot(N1, Im_rhoM_1, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,0].plot(N1, Re_rhoM_1, color="darkred",linestyle='solid', linewidth=4, label="real")
    ax[0,0].plot(N1, Im_rhoM_1, color="teal",linestyle='solid', linewidth=4, label="imag")
    
    #ax[0,0].plot(N2, Re_rhoM_2, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    #ax[0,0].plot(N2, Im_rhoM_2, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[0,0].set_ylabel("$|\\rho_{\\varphi b} \\mathcal{M}_{\\varphi b}|$")

    ax[1,0].set_ylim([10,20])
    #ax[1,0].plot(N1, diff1, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].plot(N1, diff1, color="darkorange",linestyle='solid', linewidth=4)
    
    #ax[1,0].plot(N2, diff2, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].set_ylabel("$\\Delta \\rho_{\\varphi b}$")
    ax[1,0].set_xlabel("$N$")

    ax[0,1].set_title("SA Method")
    ax[0,1].set_ylim([-0.1,1.0])
    #ax[0,1].plot(N3, Re_rhoM_3, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    #ax[0,1].plot(N3, Im_rhoM_3, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,1].plot(N3, Re_rhoM_3, color="darkred", linestyle='solid', linewidth=4, label="real")
    ax[0,1].plot(N3, Im_rhoM_3, color="teal", linestyle='solid', linewidth=4, label="imag")
    
    
    #ax[0,1].plot(N4, Re_rhoM_4, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    #ax[0,1].plot(N4, Im_rhoM_4, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[1,1].set_ylim([10,40])
    #ax[1,1].plot(N3, diff3, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,1].plot(N3, diff3, color="darkorange",linestyle='solid', linewidth=4)
    
    #ax[1,1].plot(N4, diff4, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    #ax[1,1].set_ylabel("$difference$")
    ax[1,1].set_xlabel("$N$")

    #ax[0,0].set_xscale('log')
    #ax[0,0].set_yscale('log')
    #ax[0,1].set_xscale('log')
    #ax[0,1].set_yscale('log')
    #ax[1,0].set_xscale('log')
    #ax[1,0].set_yscale('log')
    #ax[1,1].set_xscale('log')
    #ax[1,1].set_yscale('log')

    ax[0,0].legend()
    ax[0,1].legend()
    #ax[0,0].legend()
    #ax[0,0].legend()
    plt.tight_layout()
    output = "difference_between_bf_and_sa_method_upto_N500.png"
    plt.savefig(output)
    plt.show() 


#N_dependence_plots_1()   
N_dependence_plots_2() 
    