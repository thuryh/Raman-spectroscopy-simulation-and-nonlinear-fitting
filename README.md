# Raman-spectroscopy-simulation-and-nonlinear-fitting
This is a fitting tool for temperature and mole faction measurements based on the non-linear least-square-fitting method written in MATLAB.

The temperature was measured by fitting the simulated N2 Raman scattering with the experimental N2 spectra.
The mole fraction of major species (C2H4, N2, and O2) were retrieved by the measured temperature and ideal gas law.

    Code:
    get_nu_all: Caluaulate all possible Raman shift of N2 molecule. We assume there are 20 vibrational energy states and 100 rotational state (It can be 
    
    csvreader.m:  Read the data from csv files that from the experiment.
    Linearfun.m:  Linear funtion for fitting the baseline of experimental spectra.
    
