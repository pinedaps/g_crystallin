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
parser.add_argument('path', metavar='Path', type=str, nargs='?',
                    help='path where data is stored')
parser.add_argument('sigma', metavar='Sigma', type=float, nargs='?',
                    help='radius of the protein in Angtroms')
args = parser.parse_args()

data_dir = args.path
plot_dir = args.path.replace('scans','plots')

# Calculation of the B2 for a hard sphere of radius 'sigma'

sigma = args.sigma
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
plt.savefig(plot_dir+"potential.png")

fig, ax = plt.subplots(nrows=1)

b2 = []
for index, file_b2 in enumerate(files_b2):
	print('Procesing file: ', file_b2 )
	with open(data_dir+file_b2, 'r') as file:
		b2.append(json.load(file)['B2']/B2_HS)
ax.plot(T, b2, ms=5, marker='o', lw=1)
#plt.ylim(-11,5)
#plt.xlim(23,60)
plt.xlabel(r"Temperature, T [$K$]")
plt.ylabel(r"$b_2^*=b_2/b_2^{HS}$")
plt.savefig(plot_dir+"B2.png")
