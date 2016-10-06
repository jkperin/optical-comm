#!/bin/bash

# Name the job in Grid Engine
#$ -N QPSK_Analog_Delay_Penalty

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

matlab -nodesktop -singleCompThread -r "QPSK_Analog_Delay_Penalty_qsub $CPR $CPRmethod $linewidthKHz"

matlab -nodesktop -singleCompThread -r "QPSK_Analog_Delay_Penalty_qsub $CPR $CPRmethod $linewidthKHz"