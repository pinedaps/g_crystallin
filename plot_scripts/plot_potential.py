import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from math import ceil
import argparse
import os
import json
import common as cmn

#################################### Loading plotting features ##################################

cmn

######################################### Acquiring path ########################################

parser = argparse.ArgumentParser(description='Configuration of data source')
parser.add_argument('-p', '--path', type=str, 
                    help='path where the computed data is stored')
parser.add_argument('-pe', '--path_exp', type=str,
                    help='path where the experimental data is stored')
parser.add_argument('-s', '--sigma', type=float,
                    help='radius of the protein in Angtroms (default= 10 A)')
args = parser.parse_args()

data_dir = args.path
plot_dir = args.path.replace('scans','plots')
exp_dir  = args.path_exp

# Calculation of the B2 for a hard sphere of radius 'sigma'

if args.sigma:
	sigma = args.sigma
else:
	sigma = 10
	print('Protein radius not provided, a default radius of 1 nm was used!')
B2_HS = 2/3*np.pi*sigma**3
print('B2_HS = ' +str(B2_HS))

# Extract T values from filename and create the list of files 

T         = []
files_pot = []
files_b2  = []
for name in sorted(os.listdir(data_dir),reverse=True):
	if name.endswith(".dat"):
		files_pot.append(name)
		T.append(float(name.replace('scan_T','').replace('.dat','')))
	else:
 		files_b2.append(name)

print('Files to be processed:')
print(files_pot)

# Extract B2 values from filename and create the list of files

b2 = []
for index, file_b2 in enumerate(files_b2):
        with open(data_dir+file_b2, 'r') as file:
                b2.append(json.load(file)['B2']/B2_HS)

# Extract B2_reduced from experimental data

T_exp, b2_red_exp, b2_min, b2_max = np.loadtxt(exp_dir+os.listdir(exp_dir)[0], unpack=True)

fig, ax = plt.subplots(nrows=1)

for index, file_pot in enumerate(files_pot):
	plt.axhline(y=0, color='black', lw=1, ls='--', alpha=0.5)
	print('Procesing file: ', file_pot )
	R, F, U, C, exp = np.loadtxt(data_dir+file_pot, unpack=True)
	ax.plot(R, U, ms=5, marker='o', lw=1, label=str(T[index]))
plt.ylim(-11,5)
plt.xlim(23,60)
plt.xlabel(r"Distance, r [$\AA$]")
plt.ylabel(r"Potential, U(r) [$k_BT$]")
plt.legend(ncol=1, title='T [K]')
plt.title("Interaction Free Energy of $\gamma$B-crystallin", pad=15)
plt.savefig(plot_dir+"potential.png")

fig, ax = plt.subplots(nrows=1)

fig, ax = plt.subplots()

yerr = np.vstack([b2_min, b2_max])
ax.errorbar(T_exp, b2_red_exp, yerr=yerr, fmt='o', ms=5, lw=1, capsize=3, label="Bucciarelli $\it{et\ al.}$ 2016")
ax.plot(T, b2, ms=5, marker='o', lw=1, label="Duello")
plt.title("Reduced second virial coefficient ${b_2}^*$ \n of $\gamma$B-crystallin",pad=15)
plt.legend()
plt.xlabel(r"Temperature, T [$K$]")
plt.ylabel(r"$b_2^*=b_2/b_2^{HS}$")
plt.savefig(plot_dir+"B2.png")
