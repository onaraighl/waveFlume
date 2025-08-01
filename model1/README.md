model1



Fits a travelling wave to the data \\eta(x,t) emanating from the video analysis.



Requires the user to run /data/calculateFreeSurface.m and for the space-time array eta\_array

to be available on the command line, together with the vector of space points xx and the vector

of time points t.



Running  fit\_profile\_data.m gives non-linear least-squares fit.



Output variables:



* xx\_out - a refined vector of space points
* tt\_out - a refined vector of time points
* eta\_final - the fitted model, defined on the grid xx\_out (space) and tt\_out (time)
* eta\_final\_orig - the fitted model, defined on the original grid xx (space) and t (time)
* h0, a1, omega,phi1,k1 - the fitted parameters, where the model is



\\eta(x,t)=h0+A1\*cos(omega t - kx + \\phi)



Running the matlab SCRIPT fit\_profile\_data\_bootstrap generates the 2.5% and 97.5% confidence intervals for

the fitted parameters.

