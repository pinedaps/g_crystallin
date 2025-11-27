import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from math import ceil
import argparse
import os
import common as cmn

#################################### Loading plotting features ##################################

cmn

######################################### Acquiring path ########################################

parser = argparse.ArgumentParser(description='Configuration of data source')
parser.add_argument('path', metavar='Path', type=str, nargs='?',
                    default='../results/dat/',
                    help='path where data is stored (default: ../results/dat/)')
parser.add_argument('plot_setup', metavar='plot_colors_markers', type=str, nargs='?',
                    default='standard',
                    help='Colors and markers for the plot (default: standar = random colors and o marker)')
args = parser.parse_args()

plt_setup = args.plot_setup

data_dir = args.path

# Create the list of files

files = sorted(np.array( [name for name in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, name))]), reverse=True)

print('Files to be processed:')
print(files)

# Create the list of colors and fillstyles to facilitate comparison 

colors  = ['black','black','black','blue','green','red','green','red','red','red','red','orange','purple', 'purple'] 
markers = ['^',    '*',    'o',    'o',   'o',    'o',  '*',    '*',  'v',  's',  '>',  'o',     '*'    , 'o']

fig, ax = plt.subplots(nrows=1)

for index,file in enumerate(files):
	plt.axhline(y=0, color='black', lw=1, ls='--', alpha=0.5)
	print('Procesing file: ', file )
	R, F, U, C, exp = np.loadtxt(data_dir+file, unpack=True)
	if plt_setup == 'standard':
		ax.plot(R, U, ms=5, marker='o', lw=1, label=file.replace('pmf_','').replace('.dat',''))
	else:	
		ax.plot(R, U, ms=5, marker=markers[index], lw=1, color=colors[index], label=file.replace('pmf_','').replace('.dat',''))
plt.ylim(-11,5)
plt.xlim(23,60)
plt.xlabel(r"Distance,  r [$\AA$]")
plt.ylabel(r"Potential,  U [$k_BT$]")
plt.legend(ncol=1, title='T [K]')
plt.savefig("./tmp_fig/potential.png")
