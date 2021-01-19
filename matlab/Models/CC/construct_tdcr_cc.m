
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga


function [ T1_cc,T2c_cc] = construct_tdcr_cc( var_cc )
    
    kappa1 = var_cc(1,1);
    phi1 = var_cc(2,1);
    l1 = var_cc(3,1);
    kappa2 = var_cc(1,2);
    phi2 = var_cc(2,2);
    l2 = var_cc(3,2);
    
    % task space
    sect_points=50; %number of points per section

    % section 1
    [T1_cc] = trans_mat_cc(kappa1, phi1, l1, sect_points);
    T1_tip=reshape(T1_cc(length(T1_cc),:),4,4);

    % section 2
    [T2] = trans_mat_cc(kappa2, phi2-phi1, l2, sect_points);
    [T2c_cc,~] = coupletransformations(T2,T1_tip);
    T2_tip=reshape(T2c_cc(length(T2c_cc),:),4,4);
end

