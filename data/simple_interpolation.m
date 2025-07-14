function [h,a,k,phi,xxi,eta,xx,yy]=simple_interpolation(filename)

T = readtable(filename);

xxi = double(T{:,1});
eta = double(T{:,2});

% Convert to SI:

xxi=xxi/100;
eta=eta/100;

% Account for measurement error:
% Length of test secton of wave tank is 81.6 cm and not 61.6 cm which is
% what I have put into the Engauge digitizer.

xxi=xxi*(81.6/61.6);

% Generate fitted curve:

% xx=-0.2:0.01:0.8;
xx=0:0.01:0.81;
yy=interp1(xxi,eta,xx,'spline');

% There is also a small error in the digitized y-coordinates that needs to
% be corrected.

yy=(0.05/0.0552)*yy;

h=0;
a=0;
k=0;
phi=0;



end