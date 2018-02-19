#!/bin/bash

#SBATCH --job-name=GN_coeff
#SBATCH --output=GN_coeff_%A_%a.out
#SBATCH --error=GN_coeff_%A_%a.err
#SBATCH --array=1-6
#SBATCH --time=1-0
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=8
#SBATCH --qos=long

######################
# Begin work section #
######################

# Print this sub-job's task ID
echo "SLURM_ARRAY_TASK_ID is " $SLURM_ARRAY_TASK_ID

# Run your code based on the SLURM_ARRAY_TASK_ID
module load matlab

srun matlab -nodesktop -singleCompThread -r "GN_model_coeff_slurm 50 $SLURM_ARRAY_TASK_ID"