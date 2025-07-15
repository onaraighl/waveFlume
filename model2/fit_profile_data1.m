function [xx_out,tt_out,eta_final,eta_final_orig,h0,A1,omega,phi1,k1,A2,phi2,k2,A3,phi3]=fit_profile_data1(xx,t,eta_data)

% Performs a non-linear least squares and fits the data to a
% multi-mokdal form:
%
% eta(x,t)=h0+A1*cos(omega*t-k1*x+phi1)+A2*cos(omega*t+phi2)*cos(k2*x)+A3*cos(omega*t+phi3);
%

% Set upper and lower limits for the search space.

h0_lwr=0.01;
h0_upr=0.1;

A1_lwr=0;
% measured max_eta, from observations:
max_ampl= 0.01402;
A1_upr=1.1*max_ampl;


% omega=140\pm 10 RPM:
omega_lwr=13.613568165555769;
omega_upr=15.7080;

% Constraints here:
phi1_lwr=0;
phi1_upr=pi/2;

% A posteriori constraints here:
k1_lwr=23.2711;
k1_upr=28.5599;

% The 2-component:

A2_lwr=0;
A2_upr=0.4*A1_upr;

phi2_lwr=0;
phi2_upr=2*pi;

k2_lwr=10;
k2_upr=40;

% The 3-component:

A3_lwr=0;
A3_upr=0.4*A1_upr;

phi3_lwr=0;
phi3_upr=2*pi;


lb=[h0_lwr,A1_lwr,omega_lwr,phi1_lwr,k1_lwr,A2_lwr,phi2_lwr,k2_lwr,A3_lwr,phi3_lwr];
ub=[h0_upr,A1_upr,omega_upr,phi1_upr,k1_upr,A2_upr,phi2_upr,k2_upr,A3_upr,phi3_upr];

%**************************************************************************
% Initial guess for the optimization:

% p0=(ub+lb)/2;

% r=rand;
% seed is 0.40181
% min is 0.71951

r=0.75774;

p0=r*lb+(1-r)*ub;

% p0=[0.05,0.0068,2*pi/0.25,pi];

%**************************************************************************
% Define the options for fmincon:

options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...   % Show detailed iteration output
    'Diagnostics', 'on',...        % Enable diagnostics before solving
    'OptimalityTolerance', 1e-10, ...  % Set tolerance here
    'algorithm','interior-point');

% options = optimoptions('fmincon', ...
%     'OptimalityTolerance', 1e-6);  % Set tolerance here

%**************************************************************************
% Call optimizer:

p=fmincon(@fun,p0,[],[],[],[],lb,ub,[],options);

% [p,fval,exitFlag,output] = simulannealbnd(@fun,p0);

Jmin=fun(p);

format long

display(strcat('seed is',num2str(r)))
display(strcat('min is',num2str(Jmin)))

%**************************************************************************
% Parameter values:

h0=p(1);
A1=p(2);
omega=p(3);
phi1=p(4);
k1=p(5);

A2=p(6);
phi2=p(7);
k2=p(8);

A3=p(9);
phi3=p(10);

%**************************************************************************
% Calculate model:

% eta_final=0*eta_data;

xx_out=0:0.005:max(xx);
tt_out=0:0.0025:max(t);

eta_final=zeros(length(tt_out),length(xx_out));

for i=1:length(tt_out)
    for j=1:length(xx_out)
    
        t_val=tt_out(i);
        x_val=xx_out(j);
        
        % I have to put a minus sign in front of x_val because the
        % coordinate system has been flipped.
        % This ensures that the travelling wave goes the right way.
        eta_final(i,j)=h0+A1*cos(omega*t_val-k1*(-x_val)+phi1)+...
            A2*cos(omega*t_val+phi2)*cos(k2*x_val)+...
            A3*cos(omega*t_val+phi3);
    end
end


eta_final_orig=zeros(length(t),length(xx));

for i=1:length(t)
    for j=1:length(xx)
    
        t_val=t(i);
        x_val=xx(j);
        
        % I have to put a minus sign in front of x_val because the
        % coordinate system has been flipped.
        % This ensures that the travelling wave goes the right way.
        eta_final_orig(i,j)=h0+A1*cos(omega*t_val-k1*(-x_val)+phi1)+...
            A2*cos(omega*t_val+phi2)*cos(k2*x_val)+...
            A3*cos(omega*t_val+phi3);
    end
end

%**************************************************************************
%**************************************************************************
% Cost function:

    function J=fun(p)
        h0_loc=p(1);
        A1_loc=p(2);
        omega_loc=p(3);
        phi1_loc=p(4);
        k1_loc=p(5);

        A2_loc=p(6);
        phi2_loc=p(7);
        k2_loc=p(8);
        
        A3_loc=p(9);
        phi3_loc=p(10);

        eta_model=0*eta_data;

        for i_loc=1:length(t)
            for j_loc=1:length(xx)

                t_val=t(i_loc);
                x_val=xx(j_loc);

                % I have to put a minus sign in front of x_val because the
                % coordinate system has been flipped.
                % This ensures that the travelling wave goes the right way.

                eta_model(i_loc,j_loc)=h0_loc+A1_loc*cos(omega_loc*t_val-k1_loc*(-x_val)+phi1_loc)+...
                        A2_loc*cos(omega_loc*t_val+phi2_loc)*cos(k2_loc*x_val)+...
                        A3_loc*cos(omega_loc*t_val+phi3_loc);
            end
        end

        

        J=0.5*sum(sum (eta_data-eta_model).^2);
    end

%**************************************************************************
%**************************************************************************

end