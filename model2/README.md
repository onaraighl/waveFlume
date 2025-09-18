# waveFlume/model2

Fits a mixture of a standing wave and a travelling wave to the data eta(x,t) emanating from the video analysis.

Requires the user to run /data/calculateFreeSurface.m and for the space-time array eta_array to be available on the command line, together with the vector of space points xx and the vector of time points t.

# fit_profile_data1

Running  fit_profile_data1.m gives non-linear least-squares fit.

Output variables:

* xx_out - a refined vector of space points
* tt_out - a refined vector of time points
* eta_final - the fitted model, defined on the grid xx\_out (space) and tt_out (time)
* eta_final_orig - the fitted model, defined on the original grid xx (space) and t (time)
* h0, a1, omega,phi1,k1, etc. - the fitted parameters, where the model is

eta(x,t)=h0+A1\*cos(omega\*t-k1\*x+phi1)+A2\*cos(omega\*t+phi2)\*cos(k2\*x)+A3\*cos(omega\*t+phi3);

# fit_profile_data_bootstrap1

Running the matlab SCRIPT fit_profile_data_bootstrap1 generates the 2.5% and 97.5% confidence intervals for the fitted parameters.

