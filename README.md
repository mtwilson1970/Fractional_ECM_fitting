# Fractional_ECM_fitting

Marcus Wilson 
23 July 2024
This repository includes Matlab programmes for fitting an equivalent circuit model (ECM) to data of voltage V(t) and current I(t) for charging and discharging of a battery. 

***IMPORTANT***
The raw data for is contained in .mat files (available from the author marcus.wilson@waikato.ac.nz on request - they are generally too large to include with this repository - e.g. U1drive3.mat) whose three columns are time (s), current (A) and voltage (V). 
I give a short example of a .mat file in this repository: Grydrive4.mat (select this by setting use_Grydrive = true in fit_CPE_CPE_R_as_func and other file options as 'false'). The Grydrive4.mat contains two and a bit days worth of data - and can be analyzed up to two days in length. That's not much but it uploads to github!
For the 2024 manuscript Wilson et al. "Early detection of Li-ion cell failure from current-voltage time-series data" there are much larger data files which
I can share with you on request.
****************

The code is contained in the files:

gendriven_analyze.m,      which calls the function:
fit_CPE_CPE_R_as_func,    which calls:
CF_and_R_fit_twice.m

To run the code simply launch gendrive_analyze with Matlab. 

Input parameters are mostly listed at the beginning of gendriven_analyze.m, but some (notably the selection of which file to read) is
in fit_CPE_CPE_R_as_func.

The code fits a R-CPE-CPE model to day-long segments of V(t) and I(t) data. 
(The length of segment can be varied with the day_stride parameter in fit_CPE_CPE_R_as_func. )

After running, the matlab workspace will include arrays giving the daily parameter fits (the variables sufficed as _fit, e.g. a2_fit for the alpha of the second CPE - each element corresponds to each day) and impedances at different frequencies. 
Additionally, multiple plots are given as output. 

As contained in this repository, the code will analyze 2 days of Grydrive4.mat, in one day intervals. 
