#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import common as cmn
from matplotlib.offsetbox import AnchoredText

# Loading plotting features 

cmn

# Optional imports for notebook mode
def in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() and "IPKernelApp" in get_ipython().config
    except:
        return False

# ------------------------------
# Potential definition
# ------------------------------

def U_LJ(r, epsilon=0.8368, sigma=1):
    """Lennard-Jones component."""
    sigma = 1
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

def epsilon_T(T,Tc=293,epsilon_c=0.8368, c=2.52e-2):
    """Analytical solution of Eq.15 from Wennerstrom & Lindman (https://doi.org/10.1016/j.molliq.2025.128169)"""
    return epsilon_c*(T*(1/Tc+c-c*np.log(T/Tc))-c*Tc)

def U_AH(r, epsilon_c=0.8368, T_c=293, depsilon_dT=2.52e-2, r_cutoff=2.5, lam=1.0, T=293):
    """
    Generalized Ashbaugh–Hatch potential with:
      - U_LJ(r) + epsilon(T) (1 - lambda) - lambda * U_LJ(r_cutoff)   at r <= 2^(1/6)
      - lambda * (U_LJ(r) - U_LJ(r_cutoff))                           at r >  2^(1/6)
      - 0                                                             at r >= r_cutoff
    where 
      - epsilon(T) = epsilon_c * (T * (1/Tc + c - c*np.log(T/Tc)) - c*Tc)  
    """

    # Calculate r at wich the LJ potential is minimum with the chosen sigma 
    sigma = 1  
    r_min = 2**(1/6) * sigma
    
    # Calculate the epsilon depending on the provided temperature
    epsilon = epsilon_T(T, T_c, epsilon_c, depsilon_dT)
    
    if epsilon <=0:
        YELLOW = "\033[93m"
        RESET = "\033[0m"
        print(f"{YELLOW}Warning!")
        print(f"Negative value of \u03B5 = {epsilon:.3f}, choose another temperature.{RESET}")
        exit()
    U_LJ_vals     = U_LJ(r, epsilon, sigma)
    U_LJ_r_cutoff = U_LJ(r_cutoff, epsilon, sigma)
	
    # Region definitions
    U = np.zeros_like(r)

    # Region 1: r <= r_min
    mask1 = r <= r_min
    U[mask1] = U_LJ_vals[mask1] + epsilon * (1 - lam) - U_LJ_r_cutoff

    # Region 2: r_min ≤ r < r_cutoff
    mask2 = (r > r_min) & (r < r_cutoff)
    U[mask2] = lam * (U_LJ_vals[mask2] - U_LJ_r_cutoff)

    # Region 3: r >= r_cutoff → U = 0 (already zero)
    
    return U, epsilon

# ------------------------------
# Plotting
# ------------------------------
def plot_potential(epsilon_c=0.8368, T_c=293, depsilon_dT=0.02522, r_cutoff=2.5, lam=1.0, T=293):
    rmax = r_cutoff*1.2   
    r    = np.linspace(0.8, rmax, 500)

    # Static default plot

    fig, ax = plt.subplots(figsize=(7,5))
    U_LJ  = U_AH(r, epsilon_c=0.8368, r_cutoff=2.5, lam=1.0)
    U_WCA = U_AH(r, epsilon_c=0.8368, r_cutoff=2.5, lam=0.0)
    ax.plot(r, U_LJ[0], color='black', linestyle='--', alpha=0.7, label="LJ")
    ax.plot(r, U_WCA[0], color='gray', linestyle='--', alpha=0.7, label="WCA")    
    # Anchored box for static default parameters
    
    default_text = (rf"$\mathbf{{Reference\ Parameters}}$" "\n"
		    rf"$\mathbf{{\lambda_{{LJ}}}}$ = 1   &   $\mathbf{{\lambda_{{WCA}}}}$ = 0" "\n"
		    rf"$\mathbf{{\epsilon_c}}$ = 0.8368""\n"
                    rf"$\mathbf{{T_c}}$ = 293.0" "\n"
                    rf"$\mathbf{{r_{{cutoff}}}}$ = 2.5")

    box_default = AnchoredText(default_text, loc="upper right", prop=dict(size=10), frameon=True)
    box_default.patch.set_boxstyle("round")
    box_default.patch.set_facecolor("white")
    box_default.patch.set_alpha(0.9)
    box_default.patch.set_edgecolor("black")
    box_default.txt._text.set_ha("center")
    ax.add_artist(box_default)
    
    # Interactive plot

    r    = np.linspace(0.8, rmax, 500)
    U, epsilon = U_AH(r, epsilon_c, T_c, depsilon_dT, r_cutoff, lam, T)
    plt.plot(r, U,color='blue', label="Interactive")
    plt.axvline(r_cutoff, color="red", linestyle="--", label="$r_{cutoff}$")
    plt.xlim(0.5,rmax)
    plt.ylim(-2,2)
    plt.xlabel(r"Distance, r [$\AA$]")
    plt.ylabel(r"Potential, U(r) [$k_BT$]")    
    plt.title("Ashbaugh–Hatch Potential")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    # Parameter text box
    
    param_text = (rf"$\mathbf{{Interactive\ Parameters}}$" "\n"
		  rf"$\mathbf{{\epsilon}}$ = {epsilon:.4f}""\n"
                  rf"$\mathbf{{T}}$ = {T:.1f} ""\n"
		  rf"$\mathbf{{\lambda}}$ = {lam:.2f}""\n"
                  rf"$\mathbf{{r_{{cutoff}}}}$ = {r_cutoff:.2f}")

    ax = plt.gca()
    anchored_box = AnchoredText(param_text,prop=dict(size=10),loc="lower right",frameon=True)
    anchored_box.patch.set_boxstyle("round")
    anchored_box.patch.set_facecolor("white")
    anchored_box.patch.set_alpha(0.9)
    anchored_box.patch.set_edgecolor("black")
    anchored_box.txt._text.set_ha("center")
    ax.add_artist(anchored_box)
    plt.savefig('./AH_potential.png')


# ------------------------------
# Notebook interactive mode
# ------------------------------
def run_notebook_mode():
    import ipywidgets as widgets
    from ipywidgets import interact
    
    def update(epsilonc, Tc, depsilondT, rcutoff, lam, T):
        plot_potential(epsilonc, Tc, depsilondT, rcutoff, lam, T)

    interact(update,
    epsilonc=widgets.FloatSlider(
        value=0.8368, min=0.08, max=1.6, step=0.0001,
        readout_format='.4f',
        description="ε_c"
    ),
    Tc=widgets.FloatSlider(
        value=293.0, min=273, max=350, step=1,
        description="T_c"
    ),
    depsilondT=widgets.FloatSlider(
        value=0.02522, min=0.005, max=0.045, step=0.001,
        readout_format='.4f',
        description="dε/dT"
    ),
    rcutoff=widgets.FloatSlider(
        value=2.5, min=0, max=5.0, step=0.01,
        description="r_cut"
    ),
    lam=widgets.FloatSlider(
        value=1.0, min=0.0, max=2.0, step=0.05,
        description="λ"
    ),
    T=widgets.FloatSlider(
        value=293.0, min=273, max=350, step=1,
        description="T"
    ),
)

# ------------------------------
# Script mode
# ------------------------------

def run_script_mode():
    parser = argparse.ArgumentParser(description="Plot Ashbaugh–Hatch potential")
    parser.add_argument("--epsilonc", type=float, default=0.8368, help="LJ epsilon from Calvados3 model ar Tc = 293")
    parser.add_argument("--Tc", type=float, default=293.0, help="Reference temperature Calvados3 model")
    parser.add_argument("--depsilondT", type=float, default=0.02522, help="Variance of epsilon with respect to T based on averages of thermodynamics proporties of amino acid sidechains")    
    parser.add_argument("--rcutoff", type=float, default=2.5, help="Cutoff radius")
    parser.add_argument("--lambda_val", type=float, default=1.0, help="lambda parameter")
    parser.add_argument("--T", type=float, default=293.0, help="Temperature of the system")

    args = parser.parse_args()

    plot_potential(
	epsilon_c   = args.epsilonc,
        T_c         = args.Tc,
	depsilon_dT = args.depsilondT,
	r_cutoff    = args.rcutoff,
	T           = args.T,
	lam         = args.lambda_val,
         )

# ------------------------------
# Main
# ------------------------------
if __name__ == "__main__":
    if in_notebook():
        run_notebook_mode()
    else:
        run_script_mode()

