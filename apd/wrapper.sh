#!/bin/bash

### AWGN simulations
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" Preq_vs_gain_qsub.sh

# ### ISI simulations
# # ka = 0.1
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh

# # ka = 0.2
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh

# # 8-PAM
# # ka = 0.1
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh

# # ka = 0.2
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" Preq_vs_gain_qsub.sh

# #### Fiber simulations
# ka = 0.1
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

# ka=0.2
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="0" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="1" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="2" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="3" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="4" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="5" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="6" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="7" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="8" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="9" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="10" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="11" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="12" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="13" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="14" Preq_vs_gain_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30",Lkm="15" Preq_vs_gain_L_qsub.sh