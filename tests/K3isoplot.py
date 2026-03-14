import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# --- Matplotlib / LaTeX style ---
plt.rcParams.update({"font.size": 14})
plt.rc("font", family="serif", serif=["Computer Modern Roman"])
plt.rc("text", usetex=True)

# --- Given parameters ---
K3iso0, K3iso1 = 110027.81552312648, -766027.694222162
K3iso0_err, K3iso1_err = 139710.89780849387, 220542.52114499235

# Masses in lattice units
a_tmK = 0.09698
a_tmpi = 0.06906

# Thresholds (KKπ to KKππ)
E_thr_KKpi = 2.0 * a_tmK + a_tmpi
E_thr_KKpipi = 2.0 * a_tmK + 2.0 * a_tmpi

# Energy grid
E = np.linspace(E_thr_KKpi, E_thr_KKpipi, 400)

# Reference mass scale M (use KKπ threshold)
M = E_thr_KKpi

# Model: K3 = K0 + ((E^2 - M^2)/M^2) * K1
A = (E**2 - M**2) / (M**2)
K3_central = K3iso0 + A * K3iso1

# Error propagation assuming uncorrelated parameters
K3_sigma = np.sqrt(K3iso0_err**2 + (A * K3iso1_err) ** 2)

# Divide y-data by 1e5
scale = 1e5
K3_central_s = K3_central / scale
K3_sigma_s = K3_sigma / scale

# LaTeX tick formatter (all tick labels in math mode)
def latex_tick(x, pos):
    s = f"{x:.3g}"
    return rf"${s}$"

formatter = FuncFormatter(latex_tick)

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(E, K3_central_s, linewidth=1.5)
ax.fill_between(
    E,
    K3_central_s - K3_sigma_s,
    K3_central_s + K3_sigma_s,
    alpha=0.25,
    linewidth=0,
)

ax.axvline(E_thr_KKpi, linestyle="--", linewidth=1)
ax.axvline(E_thr_KKpipi, linestyle=":", linewidth=1)

# Axis labels (NO title, NO legend)
ax.set_xlabel(r"$a_t E_{\mathrm{cm}}$")
ax.set_ylabel(r"$a_t^{-2}\,\mathcal{K}_{3,\mathrm{iso}}\times 10^{-5}$")

# Tick formatting
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)

fig.tight_layout()
output = "K3isoplot.pdf"
fig.savefig(output, bbox_inches="tight")
# plt.show()
