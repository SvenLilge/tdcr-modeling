
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [res] = cross_product(u,v)
%returns cross product of columns of u with columns of v
    res=[u(2,:).*v(3,:)-v(2,:).*u(3,:);u(3,:).*v(1,:)-v(3,:).*u(1,:);u(1,:).*v(2,:)-v(1,:).*u(2,:)];
end