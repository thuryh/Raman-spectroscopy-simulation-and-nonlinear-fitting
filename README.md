# Raman-spectroscopy-simulation-and-nonlinear-fitting
This is a fitting tool for temperature and mole faction measurements based on the non-linear least-square-fitting method written in MATLAB.

The temperature was measured by fitting the simulated N2 Raman scattering with the experimental N2 spectra.
The mole fraction of major species (C2H4, N2, and O2) were retrieved by the measured temperature and ideal gas law.

 # Code:
    Isolated N2 Raman scattering:
    get_nu_all: Calculate all possible Raman shift of N2 molecule. 
    newraman_spectrum.m: Calculate the intensity at all possible Raman shift.
   
    Voigt profile to simulate broadening.
    voigt.m
    fadf.m
    
    Convolution of Isolated N2 Raman scattering and Voigt profile.
    Intensity_sim.m: Calculate nomarlized simulated N2 Raman scattering. Used for temperature measurement.
    Intensity_sim_wonormalization.m: Calculate original simulated N2 Raman scattering. Used for mole fraction measurement.
    
    Temperature measurement
    csvreader.m:  Read the data from csv files that from the experiment.
    Linearfun.m:  Linear funtion for fitting the baseline of experimental spectra.
    temperaturefitting.m: Fitting to obtain the temperature.
    
    Mole fraction measurement
    csvreaderforfraction.m:  Read the data from csv files that from the experiment.
    twopoly.m:  2 Degree polynominal funtion for fitting the baseline of experimental spectra.
    fractionforvp.m: Calculate mole fraction with best-fitted temperature and ideal gas law.
   

    
