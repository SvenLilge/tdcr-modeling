
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [kappa, phi, l] = configuration(l_t, n, d)
%Mapping from tendon lengths to configuration parameters
% INPUT:
% l_t: tendon lengths of section
% n: numbers of units 
% d: distance tendon to backbone
% OUTPUT:
% kappa: curvature
% phi: direction of curvature
% l: trunk length

temp=sqrt(l_t(1)^2+l_t(2)^2+l_t(3)^2-l_t(1)*l_t(2)-l_t(1)*l_t(3)-l_t(2)*l_t(3));

if temp==0 || imag(temp)~=0
    kappa=0;
    phi=0;
    l=l_t(1);
else
    kappa=(2*temp)/(d*(l_t(1)+l_t(2)+l_t(3)));
    phi=atan2(sqrt(3)*(l_t(2)+l_t(3)-2*l_t(1)),(3*(l_t(3)-l_t(2))));
    l=((n*d*(l_t(1)+l_t(2)+l_t(3)))/temp)*asin(temp/(3*n*d));
end


end

