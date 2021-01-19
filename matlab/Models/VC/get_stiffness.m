
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [Kbtges,Kseges] = get_stiffness(k,param)
% Assign Stiffnesses based on which tubes are present
% k:number of current section

if k==1
    Kbtges=param.tube1.Kbt+param.tube2.Kbt;
    Kseges=param.tube1.Kse+param.tube2.Kse;
elseif k==2
    Kbtges=param.tube2.Kbt;
    Kseges=param.tube2.Kse;
end


end

