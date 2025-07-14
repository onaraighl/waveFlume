function [t,xx,eta_array,kt,ht]=calculateFreeSurface()

clf

%**************************************************************************
% Set up variables for reading data from files

casename = "frame";
files=dir(strcat(casename, '_*.csv'));

% Leave this as klim1=1; don't change unless you really understand what's going in.
klim1=1; 
klim2=numel(files);

% dt from the fps on the camera

dt=1/29.86; 
t=(klim1:klim2)*dt;
% nt=length(t);


%**************************************************************************

kt=0*t;
ht=0*t;

for ctr=1:klim2

    filename=files(ctr+klim1-1).name;

    % [h,a,k,phi,xxi,eta,xx,yy]=fit_sine_curve(filename);
    [h,a,k,phi,xxi,eta,xx,yy]=simple_interpolation(filename);

    eta_array(ctr,:)=yy';

    kt(ctr)=k;
    ht(ctr)=h;

end

end

