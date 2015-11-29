#!/bin/bash

# Name the job in Grid Engine
#$ -N Preq_vs_gain_L

#tell grid engine to use current directory
#$ -cwd

# Set Email Address where notifications are to be sent
#$ -M jkperin@stanford.edu

# Tell Grid Engine to notify job owner if job 'b'egins, 'e'nds, 's'uspended is 'a'borted, or 'n'o mail
#$ -m n

# Tel Grid Engine to join normal output and error output into one file 
#$ -j y

# 
module load matlab

echo "Preq_vs_gain_L_qsub $M $ka $level_spacing $BW0 $GainBW $modBW $Lkm"

# margin_vs_L_qsub(M, ka, level_spacing, BW0GHz, GainBWGHz, modBWGHz)

matlab -nodesktop -singleCompThread -r "Preq_vs_gain_L_qsub $M $ka $level_spacing $BW0 $GainBW $modBW $Lkm"
