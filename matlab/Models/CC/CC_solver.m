
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga


function [ var_cc ] = CC_solver( l_d,n,r_disk )

    %section 1
    [kappa1,phi1,l1]=configuration(l_d(1,:),n(1),r_disk);

    %section 2
    [kappa2,phi2,l2]=configuration(l_d(2,:),n(2),r_disk);
    
    var_cc = [kappa1 kappa2;phi1 phi2;l1 l2];
end

