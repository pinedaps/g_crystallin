#!/bin/bash

# Walltime HH:MM:SS
#SBATCH -t 00:45:00

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

module purge
module add foss/2022b
module load GCC/12.3.0

source ~/duello_env/bin/activate

./epsilon_analysis.sh --pH 7.1 --T 293 --emin 0.4 --emax 0.8 --estep 0.02 --pdb pdbs/1AMM

# ./epsilon_analysis.sh --pH 7.1 --T 293 --es 0.4 --pdb pdbs/1AMM

deactivate

