#!/bin/bash

# Walltime HH:MM:SS
#SBATCH -t 08:00:00

# Job name and output files

#SBATCH -J 1AMM_T_analysis
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
module add matplotlib

source ~/duello_env/bin/activate

./T_analysis.sh --pH 7.1 --epsilon 0.6281 --tmin 288 --tmax 313 --tstep 5 --pdb pdbs/1AMM_wo --output 1AMM_epsilon_0.6281_res_0.28_dr_0.1

#./T_analysis.sh --pH 7.1 --epsilon 0.6281 --temps 290,300,310 --pdb pdbs/1AMM --outdir 1AMM_epsilon_0.6281_res_0.28_dr_0.1

deactivate

