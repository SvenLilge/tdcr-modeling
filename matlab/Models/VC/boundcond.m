
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [res,F_sigma,L_sigma] = boundcond(sigma_ind,tau,F,sect,y,r,param)
%Boundary conditions due to tendon termination
% sigma_ind: index of y for tendon termination at s=sigma
% tau: tendon tension for current section
% lastsect: 1 if last section, 0 if not
% y: state variable set
% Kse,Kbt: stiffness matrices
% refconfig: undeformed reference configuration
% r: tendon locations

  %kinematic variables at section end
  v_sigma=y(sigma_ind,13:15)';
  u_sigma=y(sigma_ind,16:18)';
  R_sigma=reshape(y(sigma_ind,4:12),3,3);
  
  %forces and moments applied by tendon loads at section end
  % p_dot=R*(lie(u)*r+v)
  p1_dot_sigma=R_sigma*([-u_sigma(3)*r(2,1);u_sigma(3)*r(1,1);-u_sigma(2)*r(1,1)+u_sigma(1)*r(2,1)]+v_sigma);
  p2_dot_sigma=R_sigma*([-u_sigma(3)*r(2,2);u_sigma(3)*r(1,2);-u_sigma(2)*r(1,2)+u_sigma(1)*r(2,2)]+v_sigma);
  p3_dot_sigma=R_sigma*([-u_sigma(3)*r(2,3);u_sigma(3)*r(1,3);-u_sigma(2)*r(1,3)+u_sigma(1)*r(2,3)]+v_sigma);

  % F_sigma=F1_sigma+F2_sigma+F3_sigma;
  F_sigma=-(tau(1)/norm(p1_dot_sigma))*p1_dot_sigma...
          -(tau(2)/norm(p2_dot_sigma))*p2_dot_sigma...
          -(tau(3)/norm(p3_dot_sigma))*p3_dot_sigma;

  % L_sigma=L1_sigma+L2_sigma+L3_sigma;
  L_sigma=-tau(1)*lie(R_sigma*r(:,1))*(1/norm(p1_dot_sigma))*p1_dot_sigma...
          -tau(2)*lie(R_sigma*r(:,2))*(1/norm(p2_dot_sigma))*p2_dot_sigma...
          -tau(3)*lie(R_sigma*r(:,3))*(1/norm(p3_dot_sigma))*p3_dot_sigma;
 
  
 if sect<2 %not last section
  %kinematic variables and orientation just before and just after section end
  v_sigmamin=y(sigma_ind-1,13:15)';
  v_sigmaplu=y(sigma_ind+1,13:15)';
  u_sigmamin=y(sigma_ind-1,16:18)';
  u_sigmaplu=y(sigma_ind+1,16:18)';
  R_sigmamin=reshape(y(sigma_ind-1,4:12),3,3);
  R_sigmaplu=reshape(y(sigma_ind+1,4:12),3,3);
  
  %constitutive laws just before and just after section end: internal force, internal moment
  [Kbt,Kse] = get_stiffness(sect,param);
  n_sigmamin=R_sigmamin*Kse*(v_sigmamin-[0;0;1]);
  m_sigmamin=R_sigmamin*Kbt*u_sigmamin;
  
  [Kbt,Kse] = get_stiffness(sect+1,param);
  n_sigmaplu=R_sigmaplu*Kse*(v_sigmaplu-[0;0;1]);  
  m_sigmaplu=R_sigmaplu*Kbt*u_sigmaplu;

  %difference between internal load and external load should be minimized
  res = [n_sigmamin-n_sigmaplu-F_sigma;m_sigmamin-m_sigmaplu-L_sigma];
  
 else   %last section, sigma+ doesn't exist
  %constitutive laws at section end: internal force, internal moment
  [Kbt,Kse] = get_stiffness(sect,param);
  n_sigma=R_sigma*Kse*(v_sigma-[0;0;1]); 
  m_sigma=R_sigma*Kbt*u_sigma;
  
  %difference between internal load and external load should be minimized
  res = [n_sigma-F_sigma-F;m_sigma-L_sigma];
     
 end

end

