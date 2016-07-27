#!/bin/bash

### AWGN simulations
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" sensitivity_vs_gain_qsub.sh

# ### ISI simulations
# # ka = 0.1
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh

# # ka = 0.2
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh

# # 8-PAM
# # ka = 0.1
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh

# # ka = 0.2
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" sensitivity_vs_gain_qsub.sh

## Last simulationsqsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="20",GainBW="290",modBW="24",lamb="1270",Lkm="0" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="1" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="2" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="3" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="4" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="5" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="6" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="7" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="8" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="9" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="10" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="11" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="12" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="13" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="14" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="15" sensitivity_vs_L_qsub.sh

qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="0" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="1" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="2" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="3" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="4" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="5" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="6" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="7" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="8" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="9" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="10" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="11" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="12" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="13" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="14" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1270",Lkm="15" sensitivity_vs_L_qsub.sh

qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="0" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="1" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="2" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="3" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="4" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="5" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="6" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="7" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="8" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="9" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="10" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="11" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="12" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="13" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="14" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="equally-spaced",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="15" sensitivity_vs_L_qsub.sh

qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="0" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="1" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="2" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="3" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="4" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="5" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="6" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="7" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="8" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="9" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="10" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="11" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="12" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="13" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="14" sensitivity_vs_L_qsub.sh
qsub -v M="4",ka="0.18",level_spacing="optimized",BW0="24",GainBW="290",modBW="30",lamb="1250",Lkm="15" sensitivity_vs_L_qsub.sh


# Last simulation gain 
