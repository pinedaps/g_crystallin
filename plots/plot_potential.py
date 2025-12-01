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
parser.add_argument('plot_setup', metavar='plot_colors_markers', type=str, nargs='?',
                    default='standard',
                    help='Colors and markers for the plot (default: standar = random colors and o marker)')
args = parser.parse_args()

plt_setup = args.plot_setup

data_dir = args.path+'/dat/'
json_dir = args.path+'/json/'
plot_dir = './tmp_fig/'+args.path.split('/')[2]+'/'

# Create the list of files

# Potential

files_pot = sorted(np.array( [name for name in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, name))]), reverse=True)

# Temperatures

T = sorted(np.array( [float(name.replace('.dat','')) for name in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, name))]), reverse=True) 

# B2

files_b2 = sorted(np.array( [name for name in os.listdir(json_dir) if os.path.isfile(os.path.join(json_dir, name))]), reverse=True)

print('Files to be processed:')
print(files_pot)

# Create the list of colors and fillstyles to facilitate comparison 

colors  = ['black','black','black','blue','green','red','green','red','red','red','red','orange','purple', 'purple'] 
markers = ['^',    '*',    'o',    'o',   'o',    'o',  '*',    '*',  'v',  's',  '>',  'o',     '*'    , 'o']

fig, ax = plt.subplots(nrows=1)

for index, file_pot in enumerate(files_pot):
	plt.axhline(y=0, color='black', lw=1, ls='--', alpha=0.5)
	print('Procesing file: ', file_pot )
	R, F, U, C, exp = np.loadtxt(data_dir+file_pot, unpack=True)
	if plt_setup == 'standard':
		ax.plot(R, U, ms=5, marker='o', lw=1, label=file_pot.replace('pmf_','').replace('.dat',''))
	else:	
		ax.plot(R, U, ms=5, marker=markers[index], lw=1, color=colors[index], label=file_pot.replace('pmf_','').replace('.dat',''))
plt.ylim(-11,5)
plt.xlim(23,60)
plt.xlabel(r"Distance,  r [$\AA$]")
plt.ylabel(r"Potential,  U [$k_BT$]")
plt.legend(ncol=1, title='T [K]')
plt.savefig(plot_dir+"potential.png")

fig, ax = plt.subplots(nrows=1)

b2 = []
for index, file_b2 in enumerate(files_b2):
	print('Procesing file: ', file_b2 )
	with open(json_dir+file_b2, 'r') as file:
		b2.append(json.load(file)['B2_reduced'])
ax.plot(T, b2, ms=5, marker='o', lw=1, label=file_b2.replace('pmf_','').replace('.json',''))
#plt.ylim(-11,5)
#plt.xlim(23,60)
plt.xlabel(r"Temperature, T [$K$]")
plt.ylabel(r"$b_2^*=b_2/b_2^{HS}$")
plt.savefig(plot_dir+"B2.png")
