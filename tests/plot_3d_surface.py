
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import sys

def plot_3d_surface(input_file, output_file='3d_surface_plot.png', N=100):
    """
    Create a 3D surface plot from input file with columns: sigp sigk realM imagM

    Parameters:
    -----------
    input_file : str
        Path to input file with 4 columns (sigp, sigk, realM, imagM)
    output_file : str
        Path to save output plot (default: '3d_surface_plot.png')
    N : int
        Number of grid points for interpolation (default: 100)
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data from file
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    # Create interpolation grid
    xi = np.linspace(sigp.min(), sigp.max(), N)
    yi = np.linspace(sigk.min(), sigk.max(), N)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate realM onto grid
    Zi = scipy.interpolate.griddata((sigp, sigk), realM, (Xi, Yi), method='linear')

    # Create 3D figure
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Create surface plot
    surf = ax.plot_surface(Xi, Yi, Zi, cmap='viridis', 
                           linewidth=0, antialiased=True, alpha=0.9)

    # Labels and title
    ax.set_xlabel(r'$\sigma_p$', fontsize=16, labelpad=10)
    ax.set_ylabel(r'$\sigma_k$', fontsize=16, labelpad=10)
    ax.set_zlabel(r'Re($M$)', fontsize=16, labelpad=10)
    ax.set_title(r'3D Surface Plot: $\sigma_p$ vs $\sigma_k$ vs Re($M$)', fontsize=18, pad=20)

    # Add colorbar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label(r'Re($M$)', fontsize=14)

    # Set viewing angle
    ax.view_init(elev=25, azim=45)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D plot saved to {output_file}")
    plt.close()


def plot_3d_with_contours(input_file, output_file='3d_surface_with_contours.png', N=100):
    """
    Create a 3D surface plot with contour projections
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    # Create grid
    xi = np.linspace(sigp.min(), sigp.max(), N)
    yi = np.linspace(sigk.min(), sigk.max(), N)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate
    Zi = scipy.interpolate.griddata((sigp, sigk), realM, (Xi, Yi), method='linear')

    # Create 3D figure
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    surf = ax.plot_surface(Xi, Yi, Zi, cmap='viridis', 
                           linewidth=0, antialiased=True, alpha=0.8)

    # Add contour lines on the surface
    ax.contour(Xi, Yi, Zi, levels=15, colors='black', linewidths=0.5, 
               linestyles='solid', offset=None)

    # Add contour projection on bottom
    zmin = np.nanmin(Zi)
    ax.contour(Xi, Yi, Zi, levels=15, colors='gray', linewidths=0.5,
               linestyles='solid', offset=zmin)

    # Labels and title
    ax.set_xlabel(r'$\sigma_p$', fontsize=16, labelpad=10)
    ax.set_ylabel(r'$\sigma_k$', fontsize=16, labelpad=10)
    ax.set_zlabel(r'Re($M$)', fontsize=16, labelpad=10)
    ax.set_title(r'3D Surface with Contours', fontsize=18, pad=20)

    # Colorbar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label(r'Re($M$)', fontsize=14)

    # Viewing angle
    ax.view_init(elev=25, azim=45)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D plot with contours saved to {output_file}")
    plt.close()


def plot_3d_scatter(input_file, output_file='3d_scatter_plot.png'):
    """
    Create a 3D scatter plot (no interpolation, just raw data points)
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    # Create 3D figure
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot with color mapping
    scatter = ax.scatter(sigp, sigk, realM, c=realM, cmap='viridis', 
                        s=20, alpha=0.6, edgecolors='black', linewidth=0.5)

    # Labels and title
    ax.set_xlabel(r'$\sigma_p$', fontsize=16, labelpad=10)
    ax.set_ylabel(r'$\sigma_k$', fontsize=16, labelpad=10)
    ax.set_zlabel(r'Re($M$)', fontsize=16, labelpad=10)
    ax.set_title(r'3D Scatter Plot (Raw Data)', fontsize=18, pad=20)

    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label(r'Re($M$)', fontsize=14)

    # Viewing angle
    ax.view_init(elev=25, azim=45)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D scatter plot saved to {output_file}")
    plt.close()


def plot_3d_wireframe(input_file, output_file='3d_wireframe_plot.png', N=50):
    """
    Create a 3D wireframe plot
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    # Create grid
    xi = np.linspace(sigp.min(), sigp.max(), N)
    yi = np.linspace(sigk.min(), sigk.max(), N)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate
    Zi = scipy.interpolate.griddata((sigp, sigk), realM, (Xi, Yi), method='linear')

    # Create 3D figure
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Wireframe plot
    ax.plot_wireframe(Xi, Yi, Zi, color='black', linewidth=0.5, alpha=0.7)

    # Labels and title
    ax.set_xlabel(r'$\sigma_p$', fontsize=16, labelpad=10)
    ax.set_ylabel(r'$\sigma_k$', fontsize=16, labelpad=10)
    ax.set_zlabel(r'Re($M$)', fontsize=16, labelpad=10)
    ax.set_title(r'3D Wireframe Plot', fontsize=18, pad=20)

    # Viewing angle
    ax.view_init(elev=25, azim=45)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"3D wireframe plot saved to {output_file}")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]

        # Optional arguments
        plot_type = sys.argv[2] if len(sys.argv) > 2 else 'surface'
        output_file = sys.argv[3] if len(sys.argv) > 3 else None
        N = int(sys.argv[4]) if len(sys.argv) > 4 else 100

        if plot_type == 'surface':
            output = output_file or '3d_surface_plot.png'
            plot_3d_surface(input_file, output, N)
        elif plot_type == 'contours':
            output = output_file or '3d_surface_with_contours.png'
            plot_3d_with_contours(input_file, output, N)
        elif plot_type == 'scatter':
            output = output_file or '3d_scatter_plot.png'
            plot_3d_scatter(input_file, output)
        elif plot_type == 'wireframe':
            output = output_file or '3d_wireframe_plot.png'
            plot_3d_wireframe(input_file, output, N)
        else:
            print(f"Unknown plot type: {plot_type}")
            print("Available types: surface, contours, scatter, wireframe")
    else:
        print("USAGE: python plot_3d_surface.py <input_file> [plot_type] [output_file] [N]")
        print()
        print("Arguments:")
        print("  input_file  : file with columns sigp sigk realM imagM")
        print("  plot_type   : surface, contours, scatter, or wireframe (default: surface)")
        print("  output_file : output PNG filename (optional)")
        print("  N           : grid resolution (optional, default: 100)")
        print()
        print("Examples:")
        print("  python plot_3d_surface.py data.txt")
        print("  python plot_3d_surface.py data.txt surface my_plot.png 150")
        print("  python plot_3d_surface.py data.txt contours")
        print("  python plot_3d_surface.py data.txt scatter")
