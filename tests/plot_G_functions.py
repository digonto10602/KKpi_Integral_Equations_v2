
import numpy as np
import matplotlib.pyplot as plt
import sys

def create_G_plots(input_file, N=10, output_prefix='G_plot'):
    """
    Create N plots of G11, G12, G21 vs k for different p values

    Parameters:
    -----------
    input_file : str
        Path to input file with columns: p k G11 G12 G21
    N : int
        Number of plots to create (default: 10)
    output_prefix : str
        Prefix for output files (default: 'G_plot')
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data from file
    p, k, G11, G12, G21 = np.loadtxt(input_file, unpack=True)

    # Get unique p values
    unique_p = np.unique(p)

    # Determine which p values to plot
    if N > len(unique_p):
        print(f"Warning: Requested {N} plots but only {len(unique_p)} unique p values found.")
        print(f"Creating {len(unique_p)} plots instead.")
        N = len(unique_p)
        p_values = unique_p
    else:
        # Select N evenly spaced p values
        indices = np.linspace(0, len(unique_p)-1, N, dtype=int)
        p_values = unique_p[indices]

    print(f"Creating {N} plots for p values: {p_values}")

    # Create plots
    for i, p_val in enumerate(p_values):
        # Filter data for this p value
        mask = np.isclose(p, p_val, rtol=1e-10)
        k_subset = k[mask]
        G11_subset = G11[mask]
        G12_subset = G12[mask]
        G21_subset = G21[mask]

        # Sort by k for better line plotting
        sort_idx = np.argsort(k_subset)
        k_sorted = k_subset[sort_idx]
        G11_sorted = G11_subset[sort_idx]
        G12_sorted = G12_subset[sort_idx]
        G21_sorted = G21_subset[sort_idx]

        # Create plot with fixed size
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))

        # Plot all three G functions vs k
        ax.plot(k_sorted, G11_sorted, linewidth=3, color='darkred', label=r'$G_{11}$')
        ax.plot(k_sorted, G12_sorted, linewidth=3, color='darkblue', label=r'$G_{12}$')
        ax.plot(k_sorted, G21_sorted, linewidth=3, color='darkgreen', label=r'$G_{21}$')

        # Labels and title
        ax.set_xlabel(r'$k$', fontsize=18)
        ax.set_ylabel(r'$G_{ij}$', fontsize=18)
        ax.set_title(r'$p = {:.4f}$'.format(p_val), fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Add legend
        ax.legend(fontsize=16, frameon=True, loc='best')

        # Remove grid lines
        ax.grid(False)

        # Save plot
        output_file = f"{output_prefix}_{i+1}.pdf"
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight')
        print(f"  Saved: {output_file} (p = {p_val:.4f})")
        plt.close()

    print(f"\nAll {N} plots created successfully!")


def create_separate_G_plots(input_file, N=10, output_prefix='G'):
    """
    Create N plots with separate subplots for G11, G12, G21

    Parameters:
    -----------
    input_file : str
        Path to input file with columns: p k G11 G12 G21
    N : int
        Number of plots to create (default: 10)
    output_prefix : str
        Prefix for output files (default: 'G')
    """

    # Set up matplotlib styling
    plt.rcParams.update({'font.size': 14})
    plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    # Read data from file
    p, k, G11, G12, G21 = np.loadtxt(input_file, unpack=True)

    # Get unique p values
    unique_p = np.unique(p)

    # Determine which p values to plot
    if N > len(unique_p):
        print(f"Warning: Requested {N} plots but only {len(unique_p)} unique p values found.")
        print(f"Creating {len(unique_p)} plots instead.")
        N = len(unique_p)
        p_values = unique_p
    else:
        # Select N evenly spaced p values
        indices = np.linspace(0, len(unique_p)-1, N, dtype=int)
        p_values = unique_p[indices]

    print(f"Creating {N} plots with 3 subplots each for p values: {p_values}")

    # Create plots
    for i, p_val in enumerate(p_values):
        # Filter data for this p value
        mask = np.isclose(p, p_val, rtol=1e-10)
        k_subset = k[mask]
        G11_subset = G11[mask]
        G12_subset = G12[mask]
        G21_subset = G21[mask]

        # Sort by k
        sort_idx = np.argsort(k_subset)
        k_sorted = k_subset[sort_idx]
        G11_sorted = G11_subset[sort_idx]
        G12_sorted = G12_subset[sort_idx]
        G21_sorted = G21_subset[sort_idx]

        # Create figure with 3 subplots
        fig, axes = plt.subplots(3, 1, figsize=(10, 12))

        # G11 subplot
        axes[0].plot(k_sorted, G11_sorted, linewidth=3, color='darkred')
        axes[0].set_ylabel(r'$G_{11}$', fontsize=18)
        axes[0].set_title(r'$p = {:.4f}$'.format(p_val), fontsize=20)
        axes[0].tick_params(axis='both', which='major', labelsize=16)
        axes[0].grid(False)

        # G12 subplot
        axes[1].plot(k_sorted, G12_sorted, linewidth=3, color='darkblue')
        axes[1].set_ylabel(r'$G_{12}$', fontsize=18)
        axes[1].tick_params(axis='both', which='major', labelsize=16)
        axes[1].grid(False)

        # G21 subplot
        axes[2].plot(k_sorted, G21_sorted, linewidth=3, color='darkgreen')
        axes[2].set_xlabel(r'$k$', fontsize=18)
        axes[2].set_ylabel(r'$G_{21}$', fontsize=18)
        axes[2].tick_params(axis='both', which='major', labelsize=16)
        axes[2].grid(False)

        # Save plot
        output_file = f"{output_prefix}_{i+1}.pdf"
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight')
        print(f"  Saved: {output_file} (p = {p_val:.4f})")
        plt.close()

    print(f"\nAll {N} plots created successfully!")


def show_available_p(input_file):
    """
    Display all available p values in the data file
    """
    p, k, G11, G12, G21 = np.loadtxt(input_file, unpack=True)
    unique_p = np.unique(p)

    print(f"\nAvailable p values ({len(unique_p)} total):")
    print("="*60)
    for i, val in enumerate(unique_p):
        # Count data points for this p
        count = np.sum(np.isclose(p, val, rtol=1e-10))
        print(f"  {i+1:3d}. p = {val:.6f}  ({count} data points)")
    print("="*60)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: python plot_G_functions.py <input_file> [mode] [N] [output_prefix]")
        print()
        print("Arguments:")
        print("  input_file    : file with columns p k G11 G12 G21")
        print("  mode          : 'combined' or 'separate' (default: 'combined')")
        print("                  or 'list' to show available p values")
        print("  N             : number of plots to create (default: 10)")
        print("  output_prefix : prefix for output files (default: 'G_plot' or 'G')")
        print()
        print("Modes:")
        print("  combined  : All G11, G12, G21 on same plot")
        print("  separate  : Three subplots stacked vertically")
        print()
        print("Examples:")
        print("  python plot_G_functions.py data.txt")
        print("  python plot_G_functions.py data.txt combined 20")
        print("  python plot_G_functions.py data.txt separate 15 my_plot")
        print("  python plot_G_functions.py data.txt list")
        sys.exit(1)

    input_file = sys.argv[1]

    # Default values
    mode = 'combined'
    N = 10
    output_prefix = None

    # Parse arguments
    if len(sys.argv) > 2:
        if sys.argv[2].lower() == 'list':
            show_available_p(input_file)
            sys.exit(0)
        mode = sys.argv[2].lower()

    if len(sys.argv) > 3:
        N = int(sys.argv[3])

    if len(sys.argv) > 4:
        output_prefix = sys.argv[4]

    # Set default output prefix based on mode
    if output_prefix is None:
        output_prefix = 'G_plot' if mode == 'combined' else 'G'

    # Create plots based on mode
    if mode == 'combined':
        create_G_plots(input_file, N, output_prefix)
    elif mode == 'separate':
        create_separate_G_plots(input_file, N, output_prefix)
    else:
        print(f"Unknown mode: {mode}")
        print("Use 'combined' or 'separate'")
        sys.exit(1)
