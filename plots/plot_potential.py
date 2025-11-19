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
args = parser.parse_args()
data_dir = args.path

# Create the list of files

files = sorted(np.array( [name for name in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, name))]), reverse=True)

# Create the list of colors and fillstyles to facilitate comparison 

colors  = ['black','black','black','blue','green','red','green','red','red','red','red','orange','purple', 'purple'] 
markers = ['^',    '*',    'o',    'o',   'o',    'o',  '*',    '*',  'v',  's',  '>',  'o',     '*'    , 'o']

fig, ax1 = plt.subplots(nrows=1)

for index,file in enumerate(files):
	print('Procesing file: ', file )
	R, F, U, C, exp = np.loadtxt(data_dir+file, unpack=True)
	ax1.plot(R, U, ms=5, marker=markers[index], fillstyle='none', lw=1, color=colors[index], label=file.replace('pmf_','').replace('.dat',''))
	
plt.ylim(-50,50)
plt.xlim(23,60)
plt.xlabel(r"Radius [$\AA$]")
plt.ylabel(r"Potential, U [$k_BT$]")
plt.legend(ncol=2)
plt.savefig("./tmp_fig/potential.png")
