
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [s,p] = construct_tdcr_ccsub(l,p_tendon,n_disk,var)

%% visualisation of the results
n_tendon = size(p_tendon,2);
p=[0;0;0];
s = 0;
for ss_i=1:n_disk
    beta=var(3*ss_i-2);
    gamma=var(3*ss_i-1);
    k=sqrt(beta^2+gamma^2);
    phi=atan2(gamma,beta);
    r=1/k;
    ni = 30;
    si=0.0:l(ss_i)/(ni-1):l(ss_i);
    p_i=r*[(1-cos(si/r))*cos(phi);(1-cos(si/r))*sin(phi);sin(si/r)]; %coordinates of points along subsegment i
    R_i=rotation_mat_ccsub(var,l,ss_i-1,0);
    p_plot=R_i*p_i+repmat(p(:,end),1,length(si)); % coordinates of end of ss_i

    if ss_i==1
        p= p_plot;
        s = si;
    else
        p=[p p_plot];
        s = [s s(end)+si];
    end

end
end