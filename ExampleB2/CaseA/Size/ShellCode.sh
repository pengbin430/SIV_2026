#!/bin/sh
#SBATCH --job-name=ShellCode
#SBATCH --partition=comp
#SBATCH --time=167:59:00
#SBATCH --mem=5000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bin.peng@monash.edu
#SBATCH --output=ShellCode.out

module load matlab/r2023a

MCR_CACHE_ROOT=$TMPDIR
export MCR_CACHE_ROOT

# Please note that MATLABROOT needs to be defined by the user in MonV2
# In MonARCH V1 this was placed in the module file
# Please point it to the appropiate directory for your version of Matlab
#
export MATLABROOT=/usr/local/matlab/r2023a

echo 10 | ./run_ShellCode.sh  $MATLABROOT