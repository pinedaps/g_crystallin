import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from math import ceil
import argparse
import os

############################# Loading the common plotting parameters #############################

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 15,
    "axes.titlesize": 17,
    "axes.labelsize": 17,
    "legend.fontsize": 13,
    "legend.edgecolor": 'none',
    "legend.frameon": True,
    "legend.framealpha": 0.5,
    "font.sans-serif": ["Helvetica"],
    "axes.facecolor": '#ffffff',
    "figure.autolayout": True,
        })

###################################################################################################

parser = argparse.ArgumentParser(description='Configuration of data source')
parser.add_argument('path', metavar='Path', type=str, nargs='?',
                    default='../results/dat/',
                    help='path where data is stored (default: ../results/dat/)')
args = parser.parse_args()

data_dir = args.path

# Create the list of files

files = np.array( [name for name in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, name))])
print(files)

exact = True
angle_res = 0.005  # radians

R, F, U, C, exp = np.loadtxt(data_dir, unpack=True)
print(R)
exit()

if exact:
    pot = pot_exact
else:
    pot = pot_interpolated

fig, ax1 = plt.subplots(nrows=1)
maxpot = max(abs(pot.min()), abs(pot.max()))

ngridx = ceil((x.max() - x.min()) / angle_res)
ngridy = ceil((y.max() - y.min()) / angle_res)
xi = np.linspace(x.min(), x.max(), ngridx)
yi = np.linspace(y.min(), y.max(), ngridy)
zi = griddata((x, y), pot, (xi[None, :], yi[:, None]), method="linear")
plt.contourf(xi, yi, zi, 20, cmap=plt.cm.RdBu)

# Overlay exact potential from each vertex (colored circles)
x, y, z, theta, phi, pot = np.loadtxt("pot_at_vertices.dat", unpack=True)
plt.scatter(
    theta, phi, c=pot, cmap=plt.cm.RdBu, s=15, marker="o", edgecolor="k", linewidths=0.1
)
plt.colorbar()  # draw colorbar
plt.clim(-maxpot, maxpot)
plt.ylim(0.0, 2.0 * np.pi)
plt.xlim(0.0, np.pi)
plt.xlabel(r"Polar angle, $\theta$")
plt.ylabel(r"Azimuthal angle, $\phi$")
plt.title(
    "Electric potential around a patchy particle:\n(Circles = ref. points on subdivided icosahedron)",
    fontsize=12,
)
plt.savefig("potential.pdf")
