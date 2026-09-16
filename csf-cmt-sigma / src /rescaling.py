"""
csf-cmt-sigma: affine rescaling test for lattice QCD symmetry energy.

Core relation (to be tested):
    sigma(T, Phi) = sigma_lat( T / (1 + beta*Phi/c^2) )

Neutron-star benchmark: beta*Phi/c^2 = 0.096  (T_c shift +9.6%)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------
BETA_PHI_C2 = 0.096          # beta*Phi/c^2 for neutron-star environment
DATA_FILE = "data/lattice_sigma_Tc.dat"
FIGURE_FILE = "figures/fig_rescaling.png"

# ------------------------------------------------------------------
# Load lattice data: two columns, T/T_c  and  sigma/sigma_0
# Accepts space-, tab-, or comma-separated plain text.
# ------------------------------------------------------------------
def load_lattice_data(path=DATA_FILE):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Data file not found: {path}\n"
            f"Expected format: two columns -> T/T_c   sigma/sigma_0"
        )
    data = np.loadtxt(path, delimiter=None, comments="#")
    T_ratio = data[:, 0]          # T / T_c(Phi=0)
    sigma_ratio = data[:, 1]      # sigma / sigma_0
    return T_ratio, sigma_ratio

# ------------------------------------------------------------------
# Affine rescaling: only T_c shifts, shape unchanged
# ------------------------------------------------------------------
def rescale_T(T_ratio, beta_phi_c2=BETA_PHI_C2):
    """Map zero-potential T to finite-Phi effective temperature."""
    return T_ratio * (1.0 + beta_phi_c2)

def rescale_sigma_curve(T_ratio, sigma_ratio, beta_phi_c2=BETA_PHI_C2):
    """
    The rescaled curve: at finite Phi the SAME sigma value is reached
    at a HIGHER T.  Equivalently, plot sigma_lat(T / (1+beta*Phi/c^2)).
    """
    return T_ratio / (1.0 + beta_phi_c2), sigma_ratio

# ------------------------------------------------------------------
# Scaling-law self-consistency check (the zero-shape-parameter claim):
#     sigma(T/T_c(Phi), Phi) = sigma_lat(T/T_c(Phi))
# This must hold EXACTLY for the rescaling ansatz.  Any deviation here
# signals a numerical/load error, not physics.
# ------------------------------------------------------------------
def scaling_law_check(T_ratio, sigma_ratio, beta_phi_c2=BETA_PHI_C2, npts=200):
    T_phi = np.linspace(T_ratio.min(), T_ratio.max(), npts)
    s_lat = np.interp(T_phi, T_ratio, sigma_ratio)               # lattice
    s_rescaled = np.interp(T_phi * (1 + beta_phi_c2),            # finite-Phi arg
                           rescale_sigma_curve(T_ratio, sigma_ratio, beta_phi_c2)[0],
                           sigma_ratio)
    max_dev = np.max(np.abs(s_lat - s_rescaled))
    return max_dev

# ------------------------------------------------------------------
# Figure: original vs rescaled
# ------------------------------------------------------------------
def make_figure(T_ratio, sigma_ratio, beta_phi_c2=BETA_PHI_C2):
    T_rescaled, sigma_rescaled = rescale_sigma_curve(T_ratio, sigma_ratio, beta_phi_c2)

    plt.figure(figsize=(6, 4.5))
    plt.plot(T_ratio, sigma_ratio, "o-", color="black",
             label=r"lattice, $\Phi=0$")
    plt.plot(T_rescaled, sigma_rescaled, "s--", color="tab:red",
             label=rf"rescaled, $\beta\Phi/c^2={beta_phi_c2}$")
    plt.xlabel(r"$T/T_c(\Phi)$")
    plt.ylabel(r"$\sigma/\sigma_0$")
    plt.legend()
    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)
    plt.savefig(FIGURE_FILE, dpi=300)
    print(f"Figure saved: {FIGURE_FILE}")

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == "__main__":
    T_ratio, sigma_ratio = load_lattice_data()

    # sanity: data should be ordered by T
    if not np.all(np.diff(T_ratio) > 0):
        raise ValueError("T/T_c column must be strictly increasing.")

    dev = scaling_law_check(T_ratio, sigma_ratio)
    print(f"[scaling-law self-check] max deviation = {dev:.3e} (must be ~0)")

    # effective critical temperature shift
    Tc_shift = 1.0 + BETA_PHI_C2
    print(f"[neutron-star env] T_c -> {Tc_shift:.4f} T_c  (+{100*BETA_PHI_C2:.1f}%)")

    make_figure(T_ratio, sigma_ratio)
