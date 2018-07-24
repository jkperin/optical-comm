#  Erbium-doped fiber amplifiers and system capacity optimization
> Developed by Jose Krause Perin during graduate school at Stanford University.

* J. Krause Perin, and J. Kahn, "Importance of Amplifier Physics in Maximizing the Capacity of Submarine Links," arXiv, 2018. [PDF](https://arxiv.org/abs/1803.07905)

## Folder description
- [edfa/](https://github.com/jkperin/optical-comm/tree/master/edfa): analysis and simulations of erbium-doped fiber amplifiers for electrical-power-limited submarine links. EDFAs are modeled using the Standard Confined-Doping (SCD) model and capacity optimization is performed using the particle swarm optimization algorithm.
  - [edfa/doc](https://github.com/jkperin/optical-comm/tree/master/edfa/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [edfa/data](https://github.com/jkperin/optical-comm/tree/master/edfa/data): data for some erbium-doped fibers
  - [edfa/f](https://github.com/jkperin/optical-comm/tree/master/edfa/f): auxiliary functions and classes used in edfa/
  - [edfa/results](https://github.com/jkperin/optical-comm/tree/master/edfa/results): folder reserved for saving files of simulations on cluster and posterior processing.
  - [edfa/validation](https://github.com/jkperin/optical-comm/tree/master/edfa/validation): scripts used to test and validate functions in edfa/

## Getting started
- Run the script [setup_edfa.m](https://github.com/jkperin/optical-comm/blob/master/edfa/setup_edfa.m). This script executes two tasks:
  - Compiles functions [GN_model_noise.cpp](https://github.com/jkperin/optical-comm/blob/master/f/GN_model_noise.cpp) and [GN_model_noise_gradient.cpp](https://github.com/jkperin/optical-comm/blob/master/f/GN_model_noise_gradient.cpp) from ../f/, which are necessary to compute the strength of the Gaussian noise model. This requires that a C++ compiler be set up in the computer. Tru,  `>> mex -setup`
  - Saves .mat files containing EDF data used by class EDF.

