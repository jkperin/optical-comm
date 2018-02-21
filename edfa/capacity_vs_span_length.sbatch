#!/bin/bash

#SBATCH --job-name=capacity_vs_span_length
#SBATCH --output=capacity_vs_span_length_%A_%a.out
#SBATCH --error=capacity_vs_span_length_%A_%a.err
#SBATCH --array=1-15
#SBATCH --time=1-0
#SBATCH --partition=normal
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

srun matlab -nodesktop -singleCompThread -r "capacity_vs_span_length_slurm 50 14350 $SLURM_ARRAY_TASK_ID"