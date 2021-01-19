
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [xmat] = lie(xvec)
%Converts a vector xvec (element of R3) to a matrix xmat (element of so(3))

xmat= [0        -xvec(3) xvec(2); ...
       xvec(3)  0       -xvec(1); ...
       -xvec(2) xvec(1) 0];

end

