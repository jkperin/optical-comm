#  Optical communications systems simulator developed by Jose Krause Perin during graduate school at Stanford University
> A brief description of the functions, scripts, and classes used in this project. More detailed comments can be found in each file.

This project contains codes for analyses and simulations of optical communications systems for data center applications. They were used for the generation of most of the results in the following publications:

1. J. K. Perin, A. Shastri, J. M. Kahn, “Design of Low-Power DSP-Free Coherent Receiver for Data Center Links,” subm. to J. of Lightwave Technol., 2016.
2. J. K. Perin, M. Sharif, J. M. Kahn, “Sensitivity Improvement in 100 Gbit/s-perWavelength Links using Semiconductor Optical Amplifiers or Avalanche Photodiodes,” to appear on J. of Lightwave Technol., 2016.
3. J. K. Perin, M. Sharif, J. M. Kahn, “Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Multicarrier,” J. of Lightwave Technol., vol.33, no.24, pp.5122-5132, Dec. 15, 2015.
4. M. Sharif, J. K. Perin, and J. M. Kahn, “Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Single-Carrier,” J. of Lightwave Technol., vol.33, no.20, pp.4268-4277, Oct. 15, 2015

## Folders description
- apd/: analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using avalanche photodiodes  
-- apd/doc/: documentation. Latex file containing some of the theoretical derivations and analyses
-- apd/f/: auxiliary functions used in apd/
-- apd/results/: folder reserved for saving files of simulations on cluster and posterior processing.

- coherent/: analysis and simulations of coherent and differentially coherent receivers. This includes DSP-based systems as well as systems based on analog signal processing
-- coherent/analog/: functions and classes for modeling building blocks in the coherent receiver based on analog signal process
-- coherent/analysis/: analysis scripts. These are typically oversimplified simulations to better understand some concepts
-- coherent/doc/: documentation. Latex file containing some of the results and analysis
-- coherent/DSP/: functions and classes for modeling building blocks in the coherent receiver based on digital signal process
-- coherent/f/: auxiliary functions and classes used in coherent/
-- coherent/results/: folder reserved for saving files of simulations on cluster and posterior processing.
-- coherent/validate/: scripts reserved for validation

- f/: folder of auxiliary functions and classes 

- gui/: graphic user interface developed during the internship at Maxim Integrated. No longer updated

- mpam/: class M-PAM and auxiliary functions for level spacing optimization.

- ofdm/: analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using OFDM

- figs/: folder containig relevant figures and scripts to generate figures.

- soa/: analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using semiconductor optical amplifiers. Class "SOA" is no longer used, and has been replaced by class "OpticalAmplifier" in f/.
