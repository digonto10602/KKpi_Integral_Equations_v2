
import numpy as np
import matplotlib.pyplot as plt
import sys

def create_multiple_plots(input_file, N=10, output_prefix='plot'):
    """
    Create N plots of realM vs sigk for different sigp values

    Parameters:
    -----------
    input_file : str
        Path to input file with columns: sigp sigk realM imagM
    N : int
        Number of plots to create (default: 10)
    output_prefix : str
        Prefix for output files (default: 'plot')
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data from file
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    # Get unique sigp values
    unique_sigp = np.unique(sigp)

    # Determine which sigp values to plot
    if N > len(unique_sigp):
        print(f"Warning: Requested {N} plots but only {len(unique_sigp)} unique sigp values found.")
        print(f"Creating {len(unique_sigp)} plots instead.")
        N = len(unique_sigp)
        sigp_values = unique_sigp
    else:
        # Select N evenly spaced sigp values
        indices = np.linspace(0, len(unique_sigp)-1, N, dtype=int)
        sigp_values = unique_sigp[indices]

    print(f"Creating {N} plots for sigp values: {sigp_values}")

    # Create plots
    for i, sigp_val in enumerate(sigp_values):
        # Filter data for this sigp value
        mask = np.isclose(sigp, sigp_val, rtol=1e-10)
        sigk_subset = sigk[mask]
        realM_subset = realM[mask]

        # Sort by sigk for better line plotting
        sort_idx = np.argsort(sigk_subset)
        sigk_sorted = sigk_subset[sort_idx]
        realM_sorted = realM_subset[sort_idx]

        # Create plot with fixed size
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))

        # Plot realM vs sigk - line only with darkred color
        ax.plot(sigk_sorted, realM_sorted, linewidth=3, color='darkred')

        # Labels and title
        ax.set_xlim(4,8.2)
        ax.set_xlabel(r'$\sigma_k/m_\pi^2$', fontsize=18)
        ax.set_ylabel(r'$m_\pi^2|\mathcal{M}_3^{(u,u)}|$', fontsize=18)
        ax.set_title(r'$\sigma_p/m_\pi^2 = {:.4f}$'.format(sigp_val), fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Remove grid lines
        ax.grid(False)

        # Save plot
        output_file = f"{output_prefix}_{i+1}.pdf"
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight')
        print(f"  Saved: {output_file} (sigp = {sigp_val:.4f})")
        plt.close()

    print(f"\nAll {N} plots created successfully!")


def create_plots_specific_sigp(input_file, sigp_values, output_prefix='plot'):
    """
    Create plots for specific sigp values

    Parameters:
    -----------
    input_file : str
        Path to input file
    sigp_values : list
        List of specific sigp values to plot
    output_prefix : str
        Prefix for output files
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)

    print(f"Creating {len(sigp_values)} plots for specified sigp values")

    for i, sigp_val in enumerate(sigp_values):
        # Find closest sigp value in data
        unique_sigp = np.unique(sigp)
        closest_idx = np.argmin(np.abs(unique_sigp - sigp_val))
        actual_sigp = unique_sigp[closest_idx]

        # Filter data
        mask = np.isclose(sigp, actual_sigp, rtol=1e-10)
        sigk_subset = sigk[mask]
        realM_subset = realM[mask]

        if len(sigk_subset) == 0:
            print(f"  Warning: No data found for sigp = {sigp_val}")
            continue

        # Sort by sigk
        sort_idx = np.argsort(sigk_subset)
        sigk_sorted = sigk_subset[sort_idx]
        realM_sorted = realM_subset[sort_idx]

        # Create plot with fixed size
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))

        # Plot realM vs sigk - line only with darkred color
        ax.plot(sigk_sorted, realM_sorted, linewidth=3, color='darkred')

        ax.set_xlabel(r'$\sigma_k$', fontsize=18)
        ax.set_ylabel(r'$|\mathcal{M}_3^{(u,u)}|$', fontsize=18)
        ax.set_title(r'$\sigma_p = {:.4f}$'.format(actual_sigp), fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Remove grid lines
        ax.grid(False)

        output_file = f"{output_prefix}_{i+1}.pdf"
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight')
        print(f"  Saved: {output_file} (sigp = {actual_sigp:.4f})")
        plt.close()

    print(f"\nAll plots created successfully!")


def show_available_sigp(input_file):
    """
    Display all available sigp values in the data file
    """
    sigp, sigk, realM, imagM = np.loadtxt(input_file, unpack=True)
    unique_sigp = np.unique(sigp)

    print(f"\nAvailable sigp values ({len(unique_sigp)} total):")
    print("="*60)
    for i, val in enumerate(unique_sigp):
        # Count data points for this sigp
        count = np.sum(np.isclose(sigp, val, rtol=1e-10))
        print(f"  {i+1:3d}. sigp = {val:.6f}  ({count} data points)")
    print("="*60)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: python create_realM_plots.py <input_file> [N or 'list'] [output_prefix]")
        print()
        print("Arguments:")
        print("  input_file    : file with columns sigp sigk realM imagM")
        print("  N             : number of plots to create (default: 10)")
        print("                  or 'list' to show available sigp values")
        print("  output_prefix : prefix for output files (default: 'plot')")
        print()
        print("Examples:")
        print("  python create_realM_plots.py data.txt")
        print("  python create_realM_plots.py data.txt 20")
        print("  python create_realM_plots.py data.txt 10 my_plot")
        print("  python create_realM_plots.py data.txt list")
        sys.exit(1)

    input_file = sys.argv[1]

    # Check if user wants to list available sigp values
    if len(sys.argv) > 2 and sys.argv[2].lower() == 'list':
        show_available_sigp(input_file)
        sys.exit(0)

    N = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    output_prefix = sys.argv[3] if len(sys.argv) > 3 else 'plot'

    create_multiple_plots(input_file, N, output_prefix)
