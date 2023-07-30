# 16-QAM Simulation Project

This project contains MATLAB code to simulate a four carrier 16-QAM communication system with a channel modeled by additive white Gaussian noise (AWGN).


## Files

NOTE: Some of the MATLAB (.m) files depend on each other. For proper operation, they should all be placed in the same folder on your local machine.

* **Project_Report.pdf** - Detailed report of the project.

* **main_system_simulation.m** - Script used for visualizing waveforms.

* **main_system_simulation_hardcorded_bits.m** - Same as main_system_simulation.m, but contains a line of code to hardcode bits.

* **ber_calculation_simulation.m** - Script to generate plot of bit error rates (BER) over a defined range of signal-to-noise ratios.

* **psd_averaging_simulation.m** - Script to generate power spectral density (PSD) plot averaged over a defined number of interations.

* **randbin.m** - Function that generates a vector of integers, each pseudo-randombly assigned to one of two logic values.

* **rrcpulse.m** - Function that generates samples of a root-raised-cosine pulse.

* **rrcpulse_visualizer.m** - Script for plotting raised root cosine pulses at different rolloff factors.