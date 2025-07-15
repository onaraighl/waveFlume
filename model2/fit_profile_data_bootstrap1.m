% function [xx_out,tt_out,eta_final,h0,A1,omega,phi1,k1,A2,phi2,k2,A3,phi3]=fit_profile_data1_bootstrap(xx,t,eta_data)

[t,xx,eta_array,kt,ht]=calculateFreeSurface();
eta_data=eta_array;

% *************************************************************************

[xx_out,tt_out,eta_final,eta_final_orig,h0,A1,omega,phi1,k1,A2,phi2,k2,A3,phi3]=fit_profile_data1(xx,t,eta_data);

p_best=zeros(10,1);

p_best(1)=h0;
p_best(2)=A1;
p_best(3)=omega;
p_best(4)=phi1;
p_best(5)=k1;

p_best(6)=A2;
p_best(7)=phi2;
p_best(8)=k2;

p_best(9)=A3;
p_best(10)=phi3;

% *************************************************************************

n_bootstrap=1000;

rng('default') % for reproducibility

residuals=eta_data-eta_final_orig;
N=length(t)*length(xx);
residuals=reshape(residuals,N,1);

U=zeros(n_bootstrap,10);

% *************************************************************************


for i=1:n_bootstrap
    residuals_synthetic0=randsample(residuals,N,'true');
    residuals_synthetic=reshape(residuals_synthetic0,length(t),length(xx));
    eta_data_synthetic=residuals_synthetic+eta_final_orig;
     [~,~,~,~,h0_temp,A1_temp,omega_temp,phi1_temp,k1_temp,A2_temp,phi2_temp,k2_temp,A3_temp,phi3_temp]=fit_profile_data1(xx,t,eta_data_synthetic);

    p_temp(1)=h0_temp;
    p_temp(2)=A1_temp;
    p_temp(3)=omega_temp;
    p_temp(4)=phi1_temp;
    p_temp(5)=k1_temp;

    p_temp(6)=A2_temp;
    p_temp(7)=phi2_temp;
    p_temp(8)=k2_temp;
    
    p_temp(9)=A3_temp;
    p_temp(10)=phi3_temp;

    U(i,:)=p_temp;
    display(i);
end

bootCI=prctile(U,[2.5,97.5]);

% % % p_best1=0*p_best;
% % % 
% % % for i=1:10
% % %     p_best1(i)=sum(U(:,i))/n_bootstrap;
% % % end
% % % 
% % % p=p_best1;
% % % 
% % % h0=p(1);
% % % A1=p(2);
% % % omega=p(3);
% % % phi1=p(4);
% % % k1=p(5);
% % % 
% % % A2=p(6);
% % % phi2=p(7);
% % % k2=p(8);
% % % 
% % % A3=p(9);
% % % phi3=p(10);
% % % 
% % % %**************************************************************************
% % % % Calculate model:
% % % 
% % % % eta_final=0*eta_data;
% % % 
% % % eta_final1=zeros(length(tt_out),length(xx_out));
% % % 
% % % for i=1:length(tt_out)
% % %     for j=1:length(xx_out)
% % % 
% % %         t_val=tt_out(i);
% % %         x_val=xx_out(j);
% % % 
% % %         % I have to put a minus sign in front of x_val because the
% % %         % coordinate system has been flipped.
% % %         % This ensures that the travelling wave goes the right way.
% % %         eta_final1(i,j)=h0+A1*cos(omega*t_val-k1*(-x_val)+phi1)+...
% % %             A2*cos(omega*t_val+phi2)*cos(k2*x_val)+...
% % %             A3*cos(omega*t_val+phi3);
% % %     end
% % % end