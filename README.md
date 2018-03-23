#  Research on optical communication systems
> Developed by Jose Krause Perin during graduate school at Stanford University.

This project contains code for analyses and simulations of optical communications systems. They were used for the generation of most of the results in the following publications:

* J. Krause Perin, and J. Kahn, "Importance of Amplifier Physics in Maximizing the Capacity of Submarine Links," arXiv, 2018. [PDF](https://arxiv.org/abs/1803.07905)
* J. Krause Perin, A. Shastri, J. Kahn, "DSP-Free Coherent Receivers for Data Center Links," OFC, 2018. [PDF](http://www.stanford.edu/~jkperin/OFC_DSP_free_coherent.pdf) 
* J. Krause Perin, A. Shastri, J. Kahn, "Data Center Links Beyond 100 Git/s per Wavelength," Photonics West, 2018. [PDF](http://www.stanford.edu/~jkperin/PW_DC_review.pdf) 
* J. Krause Perin, A. Shastri, J. Kahn, "Data Center Links Beyond 100 Git/s per Wavelength," Optical Fiber Technology, 2017. [PDF](http://www.stanford.edu/~jkperin/data_center_review.pdf) 
* J. Krause Perin, A. Shastri, and J. Kahn, "Design of Low-Power DSP-Free dsp_free_coherent Receivers for Data Center Links," J. Lightw. Technol., vol. 35, no. 21, pp. 4650–4662, 2017. [PDF](http://www.stanford.edu/~jkperin/DSP-free_coherent.pdf)
* J. Krause Perin, M. Sharif, J.M. Kahn, "Sensitivity Improvement in 100 Gbit/s-per- Wavelength Links using Semiconductor Optical Amplifiers or Avalanche Photodiodes," J. Lightw. Technol., vol. 34, no. 33, pp. 5542–5553, 2016. [PDF](http://www.stanford.edu/~jkperin/SOA_vs_APD_100G.pdf)
* J. Krause Perin, M. Sharif, J. M. Kahn, "Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Multicarrier," J. of Lightwave Technol., vol.33, no. 24, pp.5122-5132, Dec. 15, 2015. [PDF](http://ee.stanford.edu/~jmk/pubs/100.G.single-laser.multicarrier.JLT.15.pdf)
* M. Sharif, J. Krause Perin, and J. M. Kahn, "Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Single-Carrier," J. of Lightwave Technol., vol.33, no.20, pp.4268-4277, Oct. 15, 2015. [PDF](http://ee.stanford.edu/~jmk/pubs/100.G.single-laser.single-carrier.JLT.15.pdf)

## Folders description
- [edfa/](https://github.com/jkperin/optical-comm/tree/master/edfa): analysis and simulations of erbium-doped fiber amplifiers for electrical-power-limited submarine links. EDFAs are modeled using the Standard Confined-Doping (SCD) model and capacity optimization is performed using the particle swarm optimization algorithm.
  - [edfa/doc](https://github.com/jkperin/optical-comm/tree/master/edfa/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [edfa/data](https://github.com/jkperin/optical-comm/tree/master/edfa/data): data for some erbium-doped fibers
  - [edfa/f](https://github.com/jkperin/optical-comm/tree/master/edfa/f): auxiliary functions and classes used in edfa/
  - [edfa/results](https://github.com/jkperin/optical-comm/tree/master/edfa/results): folder reserved for saving files of simulations on cluster and posterior processing.
  - [edfa/validation](https://github.com/jkperin/optical-comm/tree/master/edfa/validation): scripts used to test and validate functions in edfa/

- [coherent/](https://github.com/jkperin/optical-comm/tree/master/coherent): analysis and simulations of coherent and differentially coherent receivers. This includes DSP-based systems as well as systems based on analog signal processing
  - [coherent/analog](https://github.com/jkperin/optical-comm/tree/master/coherent/analog): functions and classes for modeling building blocks in the coherent receiver based on analog signal process
  - [coherent/analysis](https://github.com/jkperin/optical-comm/tree/master/coherent/analysis): analysis scripts. These are typically oversimplified simulations to better understand some concepts 
  - [coherent/doc](https://github.com/jkperin/optical-comm/tree/master/coherent/doc): documentation. Latex file containing some of the results and analysis
  - [coherent/DSP](https://github.com/jkperin/optical-comm/tree/master/coherent/DSP): functions and classes for modeling building blocks in the coherent receiver based on digital signal process
  - [coherent/f](https://github.com/jkperin/optical-comm/tree/master/coherent/f): auxiliary functions and classes used in coherent/
  - [coherent/results](https://github.com/jkperin/optical-comm/tree/master/coherent/results): folder reserved for saving files of simulations on cluster and posterior processing.
  - [coherent/validate](https://github.com/jkperin/optical-comm/tree/master/coherent/validate): scripts reserved for validation

- [mpam/](https://github.com/jkperin/optical-comm/tree/master/mpam): M-PAM system simulation. Includes class M-PAM

- [ofdm/](https://github.com/jkperin/optical-comm/tree/master/ofdm): OFDM system simulation. Includes class OFDM, which supports DC-OFDM, ACO-OFDM, and SSB-OFDM.

- [apd/](https://github.com/jkperin/optical-comm/tree/master/apd): analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using avalanche photodiodes  
  - [apd/doc](https://github.com/jkperin/optical-comm/tree/master/apd/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [apd/f](https://github.com/jkperin/optical-comm/tree/master/apd/f): auxiliary functions used in apd/
  - [apd/results](https://github.com/jkperin/optical-comm/tree/master/apd/results): folder reserved for saving files of simulations on cluster and posterior processing.

- [f/](https://github.com/jkperin/optical-comm/tree/master/f): folder of auxiliary functions and classes.

- [figs/](https://github.com/jkperin/optical-comm/tree/master/figs): folder containig relevant figures and scripts to generate figures.

- [docs/](https://github.com/jkperin/optical-comm/tree/master/docs): foler containing documentation and github.io files.

- [stokes/](https://github.com/jkperin/optical-comm/tree/master/stokes): scripts and functions for simulation of Stokes vector receivers.

- [soa/](https://github.com/jkperin/optical-comm/tree/master/soa): analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using semiconductor optical amplifiers. Class `SOA` is no longer used and has been replaced by class `OpticalAmplifier` in [f/](https://github.com/jkperin/optical-comm/tree/master/f).