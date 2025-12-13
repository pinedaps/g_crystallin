#!/bin/bash

# Walltime HH:MM:SS
#SBATCH -t 00:05:00

# Job name and output files

#SBATCH -J pdb_T_temp
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Resources request

#SBATCH -N 1
#SBATCH --ntasks-per-node=48

# Job notification (Uncomment these line if need)
##SBATCH --mail-user=sebastian.pineda_pineda@chem.lu.se
##SBATCH --mail-type=END # other optinal types:BEGIN,END,FAIL,REQUEUE,ALL

module add foss/2022b
module add GCCcore/12.2.0
module add tbb/2021.10.0

source ~/duello_env/bin/activate

mpirun -np 48 ./T_analysis.sh --pH 7.1 --temps 290 --pdb pdbs/1AMM

deactivate

