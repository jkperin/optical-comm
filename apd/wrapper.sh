#!/bin/bash

### AWGN simulations
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="Inf",GainBW="Inf",modBW="Inf" margin_vs_gain_qsub.sh

### ISI simulations
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

# ka = 0.2
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

# ka = 0.5
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="4",ka="0.5",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

# 8-PAM
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

# ka = 0.2
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

# ka = 0.5
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_gain_qsub.sh

# qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh
# qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_gain_qsub.sh

qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh
qsub -v M="8",ka="0.5",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_gain_qsub.sh

#### Fiber simulations
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh

# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh

# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.1",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh

# ka = 0.2
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="100",modBW="30" margin_vs_L_qsub.sh

# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="200",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="20",GainBW="300",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="100",modBW="30" margin_vs_L_qsub.sh

# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh
# qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="200",modBW="30" margin_vs_L_qsub.sh

qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="equally-spaced",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh
qsub -v M="4",ka="0.2",level_spacing="optimized",BW0="30",GainBW="300",modBW="30" margin_vs_L_qsub.sh