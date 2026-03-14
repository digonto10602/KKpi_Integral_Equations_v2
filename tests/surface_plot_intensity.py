
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

def surface_plot_intensity(input_file, output_file='intensity_plot.png', N=250):
    """
    Create a surface plot from input file with columns: sig_piK sig_KK intensity

    Parameters:
    -----------
    input_file : str
        Path to input file with 3 columns (sig_piK, sig_KK, intensity)
    output_file : str
        Path to save output plot (default: 'intensity_plot.png')
    N : int
        Number of grid points for interpolation (default: 250)
    """

    # Set up matplotlib styling (matching kernel_surface_plot.py)
    plt.rcParams.update({'font.size': 18})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data from file
    # Assuming format: sig_piK  sig_KK  intensity
    sig_piK, sig_KK, intensity, garb = np.loadtxt(input_file, unpack=True)

    # Create interpolation grid
    xi = np.linspace(sig_piK.min(), sig_piK.max(), N)
    yi = np.linspace(sig_KK.min(), sig_KK.max(), N)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate data onto grid
    Zi = scipy.interpolate.griddata((sig_piK, sig_KK), intensity, (Xi, Yi), method='linear')

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))

    # Set title and labels
    ax.set_title(r'Intensity Surface Plot', fontsize=20)
    ax.set_xlabel(r'$\sigma_{\pi K}$', fontsize=20)
    ax.set_ylabel(r'$\sigma_{KK}$', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)

    # Determine contour levels based on data range
    vmin = np.nanmin(Zi)
    vmax = 1000#np.nanmax(Zi)

    # Create contour lines
    ax.contour(Xi, Yi, Zi, levels=np.linspace(vmin, vmax, num=20), 
               linewidths=0.5, colors='black', linestyles='solid')

    # Create filled contours
    h0 = ax.contourf(Xi, Yi, Zi, levels=np.linspace(vmin, vmax, num=30), cmap='RdBu_r')

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    cax.tick_params(labelsize=18)
    cbar.set_label('Intensity', fontsize=18)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


def surface_plot_with_custom_levels(input_file, output_file='intensity_plot.png', 
                                     N=250, vmin=None, vmax=None):
    """
    Create surface plot with custom intensity range

    Parameters:
    -----------
    input_file : str
        Path to input file
    output_file : str
        Output filename
    N : int
        Grid resolution
    vmin, vmax : float
        Manual intensity range (optional)
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 18})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data
    sig_piK, sig_KK, intensity, garb = np.loadtxt(input_file, unpack=True)

    # Create grid
    xi = np.linspace(sig_piK.min(), sig_piK.max(), N)
    yi = np.linspace(sig_KK.min(), sig_KK.max(), N)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate
    Zi = scipy.interpolate.griddata((sig_piK, sig_KK), intensity, (Xi, Yi), method='linear')

    # Auto-determine range if not provided
    if vmin is None:
        vmin = np.nanmin(Zi)
    if vmax is None:
        vmax = np.nanmax(Zi)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))

    # Titles and labels
    ax.set_title(r'Intensity Surface Plot', fontsize=20)
    ax.set_xlabel(r'$\sigma_{\pi K}$', fontsize=20)
    ax.set_ylabel(r'$\sigma_{KK}$', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)

    # Contours
    ax.contour(Xi, Yi, Zi, levels=np.linspace(vmin, vmax, num=20),
               linewidths=0.5, colors='black', linestyles='solid')
    h0 = ax.contourf(Xi, Yi, Zi, levels=np.linspace(vmin, vmax, num=30), cmap='viridis')

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.tick_params(labelsize=18)
    cbar.set_label('Intensity', fontsize=18)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]

        # Optional arguments
        output_file = sys.argv[2] if len(sys.argv) > 2 else 'intensity_plot.png'
        N = int(sys.argv[3]) if len(sys.argv) > 3 else 250

        surface_plot_intensity(input_file, output_file, N)
    else:
        print("USAGE: python surface_plot_intensity.py <input_file> [output_file] [N]")
        print("       input_file: file with columns sig_piK sig_KK intensity")
        print("       output_file: output PNG filename (optional, default: intensity_plot.png)")
        print("       N: grid resolution (optional, default: 250)")
