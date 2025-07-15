% function [xx_out,tt_out,eta_final,h0,A1,omega,phi1,k1,A2,phi2,k2,A3,phi3]=fit_profile_data1_bootstrap(xx,t,eta_data)

%**************************************************************************
% Seed the random-number generator for reproducibility:
rng('default') 

%**************************************************************************
% Obtain the data from fit_profile_data.m
eta_data=eta_array;

%**************************************************************************
% Obtain first estimate of model parameters:

[xx_out,tt_out,eta_final,eta_final_orig,h0,A1,omega,phi1,k1]=fit_profile_data(xx,t,eta_data);

%**************************************************************************
% Obtain first estimate of residuals:

residuals=eta_data-eta_final_orig;
N=length(t)*length(xx);
residuals=reshape(residuals,N,1);

%**************************************************************************
n_bootstrap=500;

U=zeros(n_bootstrap,5);

for i=1:n_bootstrap
    residuals_synthetic0=randsample(residuals,N,'true');
    residuals_synthetic=reshape(residuals_synthetic0,length(t),length(xx));
    eta_data_synthetic=residuals_synthetic+eta_final_orig;
    [~,~,~,~,h0_temp,A1_temp,omega_temp,phi1_temp,k1_temp]=fit_profile_data(xx,t,eta_data_synthetic);

    p_temp(1)=h0_temp;
    p_temp(2)=A1_temp;
    p_temp(3)=omega_temp;
    p_temp(4)=phi1_temp;
    p_temp(5)=k1_temp;

    U(i,:)=p_temp;
    display(i);
end

bootCI=prctile(U,[2.5,97.5]);

%**************************************************************************

