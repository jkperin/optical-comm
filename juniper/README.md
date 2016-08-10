#  Description of codes developed during internship at Juniper Networks
> A brief description of the functions, scripts, and classes used in this project. More detailed comments can be found in each file.

This project contains codes for analysis, simulation, and experiments of 4-PAM, Duobinary 4-PAM, and OFDM/DMT for optical systems targeting data center applications up 80 km. 

Many of the scripts/functions/classes in juniper/ depend on other functions/classes in the folders f/, mpam/, soa/, apd/, ofdm/, and ofdm/f/. Make sure that these folders are added to the path before running scripts.

## juniper/ folders description
- analysis: contains scripts to test ideas. Generally this simulations are overly simplified in order to better understand concepts.
- data: data from experiments. This is not pushed to the repository.
- figs: scripts use to generate figures.
- f: auxiliary functions

## BER simulation description
- BER simulations for each modulation scheme follows the same strucure.
1. First there is a main file, where all simulation parameters are defined. e.g., [juniper/pam_main.m](https://github.com/jkperin/research/juniper/pam_main.m), [juniper/dmt_main.m](https://github.com/jkperin/research/juniper/dmt_main.m).
2. The main file calls a function, which will swipe some parameter and calculate the BER. e.g., [juniper/preamplified_sys_ber.m](https://github.com/jkperin/research/juniper/preamplified_sys_ber.m), [juniper/dmt_ber.m](https://github.com/jkperin/research/juniper/dmt_ber.m). In this project the swept parameter is the launched power.
3. At each value being swept, the function calls other functions that will calculate the BER using different methods such as Montecarlo simulation, Gaussian approximation, etc. 
-- [juniper/f/ber_preamp_sys_montecarlo.m] (https://github.com/jkperin/research/juniper/f/ber_preamp_sys_montecarlo.m) - Calculates BER from Montecarlo simulation and Gaussian approximation.
-- [juniper/f/ber_preamp_sys_montecarlo_labsetup.m] (https://github.com/jkperin/research/juniper/f/ber_preamp_sys_montecarlo_labsetup.m) - It does the same thing as the function abover, but the receiver is assumed to do coherent detection (as done in the lab), as opposed to direct detection.
> These functions also estimate BER using Gaussian approximation. That is, the conditional probabilities for each level are assumed to be Gaussian whose variance depends on the intensity level. This has been proven to be fairly accurate whenever inter-symbol interference and other effects such as quantization are not dominant.
-- [juniper/f/ber_dc_ofdm.m] (https://github.com/jkperin/research/juniper/f/ber_dc_ofdm.m) - Calculates BER from Montecarlo simulation.

## Scripts used in experiments
- [juniper/generate_dac_ofdm_waveform.m](https://github.com/jkperin/research/juniper/generate_dac_ofdm_waveform.m) - Generates OFDM waveform, save it to file, and load waveform to DAC if DAC global variable is not empty.
- [juniper/generate_dac_pam_waveform.m](https://github.com/jkperin/research/juniper/generate_dac_pam_waveform.m) - Generates PAM waveform, save it to file, and load waveform to DAC if DAC global variable is not empty.
- [juniper/process_pam_waveforms.m](https://github.com/jkperin/research/juniper/process_pam_waveforms.m) - Process PAM waveforms captured from DSO. Input waveform can be either as a matrix or as the file name to the .h5 file saved by the DSO. This function returns the BER.
- [juniper/process_ofdm_waveforms.m](https://github.com/jkperin/research/juniper/process_ofdm_waveforms.m) - Process OFDM waveforms captured from DSO. Input waveform can be either as a matrix or as the file name to the .h5 file saved by the DSO. This function returns the BER.
- [juniper/live_waveform_processing.m](https://github.com/jkperin/research/juniper/live_waveform_processing.m) - Constantly captures PAM waveforms from DSO and processes them showing BER, equalizer convergence, symbols after equalization, and estimated modulator biasing.
- [juniper/capture_PAM_BER_vs_OSNR_experiment.m](https://github.com/jkperin/research/juniper/capture_PAM_BER_vs_OSNR_experiment.m) - This function is used to capture OSNR vs BER curves. After the user specifies the current OSNR, this function captures waveforms from DSO and processes it to compute the BER. The waveforms and other variables such as BER are saved once the experiment is done.
