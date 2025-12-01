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
def U_LJ(r, epsilon=1.0, sigma=1.0):
    """Lennard-Jones component."""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

def U_AH(r, sigma=1.0, epsilon_c=0.8368, T_c=293, depsilon_dT=0.02522, r_min=2**(1/6), r_cutoff=2.5, lam=1.0, T=293):
    """
    Generalized Ashbaugh–Hatch potential with:
      - U_LJ(r) + epsilon(T) (1 - lambda) - lambda * U_LJ(r_cutoff)   at r <= r_min
      - lambda * (U_LJ(r) - U_LJ(r_cutoff))                           at r >  r_min
      - 0                                                             at r >= r_cutoff
    where 
      - epsilon(T) = e_c(default=0.8368 from Calvados3 model) / Tc * (T - d(epsilon)/dT*(T-Tc)**2)  
    """

    if r_min is None:
    	r_min = 2**(1/6)*sigma 
    else:
    	0
    if r_cutoff is None:
    	r_cutoff = 2.5*sigma
    else:
    	0
    
    epsilon = epsilon_c/T_c * (T - depsilon_dT*(T-T_c)**2)
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
def plot_potential(sigma=1.0, epsilon_c=0.8368, T_c=293, depsilon_dT=0.02522, r_min=2**(1/6), r_cutoff=2.5, lam=1.0, T=293, ymax=10):
    rmax = r_cutoff*1.2   
    r    = np.linspace(0.8, rmax, 500)
    
    # Static default plot

    fig, ax = plt.subplots(figsize=(7,5))
    r_default = np.linspace(0.8, 10, 500)
    U_default = U_AH(r_default)
    ax.plot(r_default, U_default[0], color='gray', linestyle='--', alpha=0.5, label="Default")
    
    # Anchored box for static default parameters
    
    default_text = (rf"$\mathbf{{Default\ Parameters}}$" "\n"
		    rf"$\mathbf{{\epsilon_c}}$        = 0.8368""\n"
                    rf"$\mathbf{{T_c}}$        = 293.0" "\n"
		    rf"$\mathbf{{\lambda}}$        = 1.000" "\n"
                    rf"$\mathbf{{r_{{min}}}}$   = {2**(1/6):.3f}" "\n"
                    rf"$\mathbf{{r_{{cutoff}}}}$ = 2.500")

    box_default = AnchoredText(default_text, loc="upper left", prop=dict(size=10), frameon=True)
    box_default.patch.set_boxstyle("round")
    box_default.patch.set_facecolor("white")
    box_default.patch.set_alpha(0.9)
    box_default.patch.set_edgecolor("black")
    box_default.txt._text.set_ha("center")
    ax.add_artist(box_default)
    
    # Interactive plot
    
    U, epsilon = U_AH(r, sigma, epsilon_c, T_c, depsilon_dT, r_min, r_cutoff, lam, T)
    plt.plot(r, U,color='blue', label="Interactive")
    plt.axvline(r_cutoff, color="red", linestyle="--", label="$r_{cutoff}$")
    plt.xlim(0.5*sigma,rmax)
    plt.ylim(None,ymax)
    plt.xlabel(r"Distance, r [$\AA$]")
    plt.ylabel(r"Potential, U(r) [$k_BT$]")    
    plt.title("Ashbaugh–Hatch Potential")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    # Parameter text box
    
    param_text = (rf"$\mathbf{{Interactive\ Parameters}}$" "\n"
		  rf"$\mathbf{{\epsilon}}$        = {epsilon:.3f}""\n"
                  rf"$\mathbf{{T}}$        = {T:.1f}""\n"
		  rf"$\mathbf{{\lambda}}$        = {lam:.3f}""\n"
                  rf"$\mathbf{{r_{{min}}}}$   = {r_min:.3f}""\n"
                  rf"$\mathbf{{r_{{cutoff}}}}$ = {r_cutoff:.3f}")

    ax = plt.gca()
    anchored_box = AnchoredText(param_text,prop=dict(size=10),loc="upper right",frameon=True)
    anchored_box.patch.set_boxstyle("round")
    anchored_box.patch.set_facecolor("white")
    anchored_box.patch.set_alpha(0.9)
    anchored_box.patch.set_edgecolor("black")
    anchored_box.txt._text.set_ha("center")
    ax.add_artist(anchored_box)
    plt.savefig('./tmp_fig/AH_potential.png')


# ------------------------------
# Notebook interactive mode
# ------------------------------
def run_notebook_mode():
    import ipywidgets as widgets
    from ipywidgets import interact
    
    def update(sigma, epsilonc, Tc, depsilondT, rmin, rcutoff, lam, T, ymax):
        plot_potential(sigma, epsilonc, Tc, depsilondT, rmin, rcutoff, lam, T, ymax)
    interact(update,
        sigma=widgets.FloatSlider(value=1.0, min=0.1, max=5.0, step=0.1),
	epsilonc=widgets.FloatSlider(value=0.8368, min=0.1, max=5.0, step=0.001,readout_format='.4f'),
        Tc=widgets.FloatSlider(value=293.0, min=273, max=350, step=1),
	depsilondT=widgets.FloatSlider(value=0.02522, min=0.005, max=0.045, step=0.001, readout_format='.4f'),
        rmin=widgets.FloatSlider(value=1.122, min=1.122, max=5.0, step=0.1),
        rcutoff=widgets.FloatSlider(value=2.5, min=1.122, max=5.0, step=0.01),
        lam=widgets.FloatSlider(value=1.0, min=0.0, max=2.0, step=0.05),
        T=widgets.FloatSlider(value=293.0, min=273, max=350, step=1),
	ymax=widgets.FloatSlider(value=10, min=0, max=50, step=1),)

# ------------------------------
# Script mode
# ------------------------------

def run_script_mode():
    parser = argparse.ArgumentParser(description="Plot Ashbaugh–Hatch potential")
    parser.add_argument("--sigma", type=float, default=1.0, help="LJ sigma")
    parser.add_argument("--epsilonc", type=float, default=0.8368, help="LJ epsilon from Calvados3 model ar Tc = 293")
    parser.add_argument("--Tc", type=float, default=293.0, help="Reference temperature Calvados3 model")
    parser.add_argument("--depsilondT", type=float, default=0.02522, help="Variance of epsilon with respect to T based on averages of thermodynamics proporties of amino acid sidechains")    
    parser.add_argument("--rmin", type=float, default=2.0**(1/6), help="Radius where the U(r) is minimum")
    parser.add_argument("--rcutoff", type=float, default=2.5, help="Cutoff radius")
    parser.add_argument("--lambda_val", type=float, default=1.0, help="lambda parameter")
    parser.add_argument("--T", type=float, default=293.0, help="Temperature of the system")
    parser.add_argument("--ymax", type=float, default=10, help="Maximum limit in the y-axis")

    args = parser.parse_args()

    plot_potential(
        sigma       = args.sigma,
	epsilon_c   = args.epsilonc,
        T_c         = args.Tc,
	depsilon_dT = args.depsilondT,
	r_min       = args.rmin,
	r_cutoff    = args.rcutoff,
	T           = args.T,
	lam         = args.lambda_val,
	ymax        = args.ymax,
         )

# ------------------------------
# Main
# ------------------------------
if __name__ == "__main__":
    if in_notebook():
        run_notebook_mode()
    else:
        run_script_mode()

