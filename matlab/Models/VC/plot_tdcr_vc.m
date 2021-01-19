
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [Fig] = plot_tdcr_vc(g,s,l)
%Visualization of the backbone of the continuum robot for 2 sections
  fig=gcf;
  hold on

  radius1=4;
  radius2=3;
  
  lmax=l(1)+l(2);
  tube_ind1 = find(s<=l(1));
  tube_end1=tube_ind1(end);
   
  % draw backbone
  plot3(g(1:tube_end1,13), g(1:tube_end1,14),g(1:tube_end1,15),...
        'LineWidth',radius1,'Color',[0 0 1]); 
  plot3(g(tube_end1:end,13), g(tube_end1:end,14),g(tube_end1:end,15),...
        'LineWidth',radius2,'Color',[0 1 0]); 
  
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
%     view([90,0])
    daspect([1 1 1])

end

