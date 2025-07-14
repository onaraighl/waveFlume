function [h,a,k,phi,xxi,eta,xx,yy]=fit_sine_curve(filename)

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

% There is also a small error in the digitized y-coordinates that needs to
% be corrected.

eta=(0.05/0.0552)*eta;

% p0=[0.05,0.0068,2*pi/0.25,pi];

lb=[0.01,0.001    2*pi,     0];
ub=[0.1, 0.2, 2*pi/0.1, 2*pi];

% p0=(ub+lb)/2;
p0=[0.05,0.0068,2*pi/0.25,pi];

% Define the options for fmincon:

% options = optimoptions('fmincon', ...
%     'Display', 'iter-detailed', ...   % Show detailed iteration output
%     'Diagnostics', 'on',...        % Enable diagnostics before solving
%     'OptimalityTolerance', 1e-6);  % Set tolerance here

options = optimoptions('fmincon', ...
    'OptimalityTolerance', 1e-6);  % Set tolerance here


p=fmincon(@fun,p0,[],[],[],[],lb,ub,[],options);

h=p(1);
a=p(2);
k=p(3);
phi=p(4);

% Generate fitted curve:

xx=0:0.01:0.81;
yy=h+a*sin(k*xx+phi);

%**************************************************************************
% Cost function to be minimized:

    function J=fun(p)
        h_local=p(1);
        a_local=p(2);
        k_local=p(3);
        phi_local=p(4);

        y_model=h_local+a_local*sin(k_local*xxi+phi_local);

        J=0.5*sum( (eta-y_model).^2);
    end










end