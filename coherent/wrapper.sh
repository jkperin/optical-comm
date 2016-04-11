#!/bin/bash

## Modulator="SiPhotonics", BW = 40 GHz
qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="40",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

## Modulator="SiPhotonics", BW = 30 GHz
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

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="SiPhotonics",ModBW="30",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

## Modulator="MZM"
qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="DPLL",CPRtaps="0",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="7",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="4" qpsk_ber_qsub.sh

qsub -v fiberLength="0",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="1",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="2",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="3",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="4",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="5",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="6",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="7",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="8",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="9",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh
qsub -v fiberLength="10",Modulator="MZM",ModBW="0",CPRAlgorithm="Feedforward",CPRtaps="15",linewidth="200",fOffset="0",ros="5/4",ENOB="Inf" qpsk_ber_qsub.sh