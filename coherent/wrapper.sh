#!/bin/bash
#### DPLL vs feedforward 
## linewidth = 200 kHz and frequency offset = 8 GHz
qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

## Feedforward with small complexity
qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="5",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

## Feedforward with large number of taps
qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="8",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
