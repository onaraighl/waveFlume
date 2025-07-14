Data files for analysing the video.



Link to video:



https://www.youtube.com/watch?v=ahXxWVJAG-E



Data for first 120 frames have been extracted and saved as csv files - obtainable from this directory.



To extract the data into a space-time array, run:



&nbsp;\[t,xx,eta\_array,kt,ht]=calculateFreeSurface();



Here:



* t is a vector of points in time of size Nt
* xx is a vector of points in space of size Nx
* eta\_array is a space-time array of size Nt x Nx
* The other two variables kt and ht are redundant.



This .m function is looping over all the frames.  At each frame, the corresponding .csv file is being read.

The relevant values are then interpolated on to the uniform vector xx.  There are two methods for doing this,

only one is active at a time and the other is commented out.  These are:



* \[h,a,k,phi,xxi,eta,xx,yy]=fit\_sine\_curve(filename);



Performs a nonlinear least squares fit and fits the data in the csv file called "filename" to a sinusoidal form.



* \[h,a,k,phi,xxi,eta,xx,yy]=simple\_interpolation(filename);



fits the data in the csv file called "filename" using simple interpolation and cubic splines.



Here, h is interpolated free-surface profile; for each point in the array xx there is a corresponding h-value.



Each frame in the video analysis produces a vector h, these are then read into the space-time array eta\_array.

