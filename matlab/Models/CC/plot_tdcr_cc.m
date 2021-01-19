
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [ fig ] = plot_tdcr_cc( T1_cc,T2_cc,L )

    fig = gcf;
    radius1=4;
    radius2=3;
    lmax = L(1)+L(2);

    hold on;
    plot3(T2_cc(:,13),T2_cc(:,14),T2_cc(:,15),'LineWidth',radius2,'Color',[0 1 0]);
    plot3(T1_cc(:,13),T1_cc(:,14),T1_cc(:,15),'LineWidth',radius1,'Color',[0 0 1]);


    %Ground plate
  color = [1 1 1]*0.9;
  squaresize = 0.02;
  fill3([1 1 -1 -1]*squaresize,[-1 1 1 -1]*squaresize,[0 0 0 0],color);

  axis([-lmax lmax -lmax lmax 0 lmax])
  xlabel('x (m)')
  ylabel('y (m)')
  zlabel('z (m)')
  grid on
  view([0.5 0.5 0.5])
  daspect([1 1 1])
    hold off;

end

