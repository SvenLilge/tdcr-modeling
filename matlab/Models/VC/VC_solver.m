
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [y,init,exitflag] = VC_solver(init_guess,F,Ftex,L,param)


% tube-weight-distribution
q_tube=[param.tube1.q param.tube2.q];

% direction of gravity
eg= [0 0 -1]';

% tendon locations
r=[param.tendon.r1 param.tendon.r2 param.tendon.r3];

%options for fsolve
options1 = optimset('Display','iter','TolFun',1e-15,'tolx',1e-6); 
        
tic;
[init,res,exitflag,output]=fsolve(@optim_f, init_guess(1:6), options1); 
toc;

%% solver

  function [res] = optim_f(init_guess)
  % Iterate to find an init guess so that the boundary conditions at the tip
  % of the tube are satisfied

  %solve initial value problem
  [y]=run_IVP(init_guess); 
  
  ind=length(y(:,19));
  [res2] = boundcond(ind(1),F(:,2),Ftex,2,y,r,param);
 
%   res = [res1;res2];
  res = res2;
  end

  function [y] = run_IVP(init_guess)
  %Runs initial value problem

  %Break Points (section ends)
  Break_Points=[0,L(1),sum(L(1:2))];

  %initial values  
  s_init=0;
  R_init=[1,0,0;0,1,0;0,0,1];  
  p_init=[0,0,0];
  
  y=[p_init,reshape(R_init,1,9),init_guess(1:6)',s_init];
  
  %numerical integration routine
  options=odeset('AbsTol',1e-6,'RelTol',1e-6,'InitialStep',0.005);
  
  for k=1:2
    [~,b]=ode45(@(s,y) deriv(y,k),[Break_Points(k),Break_Points(k+1)],y(end,:),options); 
    if k<2

      y = [y;b];
      R = reshape(y(end,4:12),3,3);
      v = y(end,13:15)';
      u = y(end,16:18)';
      [Kbt,Kse]=get_stiffness(k,param); 
      
      ind=length(y(:,19));
      [~,F_sigma,L_sigma] = boundcond(ind(1),F(:,1),Ftex,2,y,r,param);
      
      yinit = [b(end,1:12) (v-Kse\R'*F_sigma)' (u-Kbt\R'*L_sigma)' b(end,19)]; 
      y=[y;yinit]; 
      
    else
      y=[y;b];
    end
  end
  
  y= unique(y,'rows','stable'); 
  end

  function [dy] = deriv(y,k)
  %This function contains the differential equations describing the system
  %Overview of variables
  % u: angular rate of change of g
  % v: linear rate of change of g
  % R: orientation of backbone
  
  u = [y(16);y(17);y(18)];  
  v = [y(13);y(14);y(15)];  
  R = reshape(y(4:12),3,3);

  %stiffness matrices (Kbt:bending and torsion, Kse: shear and extension)
  [Kbt,Kse]=get_stiffness(k,param); 
        
  %external force- and momentdistribution
  % fe=param.roh*sum(A_tube(k:2))*9.81*[0;0;1];
  fe=(sum(q_tube(k:2))+param.q_disks(k))*eg;
  le=[0;0;0];
    
  % intermediate matrix and vector quantities (used in differential equations)
  [A,B,G,H,c,d] = intermedquant(u,v,R,r,Kse,Kbt,F,fe,le,k);
 
  % Differential equations 
  vu_dot = [Kse + A, G; B, Kbt+H]\[d;c];
  p_dot=R*v;
%   p_dot = R*[0;0;1];
  R_dot=R*lie(u);
  s_dot=1;
    
  dy=[p_dot;reshape(R_dot,9,1);vu_dot;s_dot];

  end



end

