
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga


function [y_total,m_states] = VCref_solver(init_guess,F,Ftex,Mtex,L,param,n)

Break_Points=[0:L(1)/n(1):L(1) L(1)+L(2)/n(2):L(2)/n(2):sum(L(1:2))];
 [Kbt,Kse]=get_stiffness(1,param); 
 
%options for fsolve
options1 = optimset('Display','iter','TolFun',1e-4); 
% tendon locations
r=[param.tendon.r1 param.tendon.r2 param.tendon.r3 param.tendon.r1 param.tendon.r2 param.tendon.r3];

y_total = [];

tic;
[init,res,exitflag,output]=fsolve(@optim_f, init_guess, options1);
toc;
 
%% solver

  function [res] = optim_f(input)
  % Iterate to find an init guess so that the boundary conditions at the tip
  % of the tube are satisfied
  
  m_states=zeros(sum(n)+1,19);
  %solve initial value problem
  res_mat=zeros(6*sum(n),1);
  s_init=0;
  R_init=[1,0,0;0,1,0;0,0,1];  
  p_init=[0,0,0];
  
  init_state_segment=[p_init,reshape(R_init,1,9),input(1:6)',s_init];
  m_states(1,:)=init_state_segment;
  
  y_total = [];
  
  for k=1:sum(n)
      
    [y]=run_IVP(init_state_segment,k); 
    y_total = [y_total;y];
    m_states(k+1,:)=y(end,:);
 
    p_c = y(end,1:3)';
    R_c = reshape(y(end,4:12),3,3);
    v_c = y(end,13:15)';
    u_c = y(end,16:18)';
    
    u_n=[0;0;0]; v_n=[0;0;1];
    if k<sum(n)
        v_n=input(6*k+1:6*k+3);
        u_n=input(6*k+4:6*k+6);
    end
    n_c=R_c*Kse*(v_c-[0;0;1]);
    m_c=R_c*Kbt*u_c;
    n_n=R_c*Kse*(v_n-[0;0;1]);
    m_n=R_c*Kbt*u_n;
    res_mat((k-1)*6+1:(k-1)*6+3)=n_c-n_n;
    res_mat((k-1)*6+4:(k-1)*6+6)=m_c-m_n;
    
    if k<sum(n)
        init_state_segment=[y(end,1:12) input(6*k+1:6*k+6)' y(end,19)];%updating initial guess for next subsegment
    end
  end
  nt = length(F)/length(n); % Number of tendons per section
  Fdisk = [repmat(F,n(1),1);[zeros(n(2),nt) repmat(F(nt+1:end),n(2),1)]];  
  for k=1:sum(n)
    
    %current frame
    y = m_states(k+1,:);
    p_c = y(1:3)';
    R_c = reshape(y(4:12),3,3);

    %previous frame
    y = m_states(k,:);
    p_p = y(1:3)';
    R_p = reshape(y(4:12),3,3);
    
    tendon_pos_cur=zeros(3,6);
    tendon_pos_prev=zeros(3,6);

    for m=1:6
        tendon_pos_cur(:,m)=R_c*r(:,m)+p_c;%current frame tendon position
        tendon_pos_prev(:,m)=R_p*r(:,m)+p_p; %previous frame tendon position
    end
    
    if k<sum(n)
        %next frame
        y = m_states(k+2,:);
        p_n = y(1:3)';
        R_n = reshape(y(4:12),3,3);
        tendon_pos_next=zeros(3,6);
        for m=1:6
            tendon_pos_next(:,m)=R_n*r(:,m)+p_n;%next frame tendon position
        end
    end
    
    zi=R_c(:,3);
    F_tendon=[0;0;0];
    M_tendon=[0;0;0];
    
    if k==sum(n)
        for m=1:6
            dir1=tendon_pos_prev(:,m)-tendon_pos_cur(:,m);
            force1=Fdisk(k,m)*dir1./norm(dir1);
            F_tendon=F_tendon+force1;
            pos=tendon_pos_cur(:,m)-p_c;
            M_tendon=M_tendon+cross(pos,force1);
        end
        res_mat((k-1)*6+1:(k-1)*6+3)=res_mat((k-1)*6+1:(k-1)*6+3)-F_tendon-Ftex;
        res_mat((k-1)*6+4:(k-1)*6+6)=res_mat((k-1)*6+4:(k-1)*6+6)-M_tendon-Mtex;
    elseif (k==n(1))
        for m=1:6
            dir1=tendon_pos_prev(:,m)-tendon_pos_cur(:,m);
            dir2=tendon_pos_next(:,m)-tendon_pos_cur(:,m);
            force1=Fdisk(k,m)*dir1./norm(dir1)+Fdisk(k+1,m)*dir2./norm(dir2);

            if m>3
                force1=force1-(zi'*force1)*zi; %removing the z component of force for the second set of tendons
            end
            pos=tendon_pos_cur(:,m)-p_c;
            F_tendon=F_tendon+force1;
            M_tendon=M_tendon+cross(pos,force1);
        end
        res_mat((k-1)*6+1:(k-1)*6+3)=res_mat((k-1)*6+1:(k-1)*6+3)-F_tendon;
        res_mat((k-1)*6+4:(k-1)*6+6)=res_mat((k-1)*6+4:(k-1)*6+6)-M_tendon;
    else
        for m=1:6
            dir1=tendon_pos_prev(:,m)-tendon_pos_cur(:,m);
            dir2=tendon_pos_next(:,m)-tendon_pos_cur(:,m);
            force1=Fdisk(k,m)*dir1./norm(dir1)+Fdisk(k,m)*dir2./norm(dir2);
            force1=force1-(zi'*force1)*zi;
            F_tendon=F_tendon+force1;
            pos=tendon_pos_cur(:,m)-p_c;
            M_tendon=M_tendon+cross(pos,force1);
        end
        res_mat((k-1)*6+1:(k-1)*6+3)=res_mat((k-1)*6+1:(k-1)*6+3)-F_tendon;
        res_mat((k-1)*6+4:(k-1)*6+6)=res_mat((k-1)*6+4:(k-1)*6+6)-M_tendon;
    end
  end
  
  res = res_mat;
  end

  function [y] = run_IVP(y_init,k)
  %Runs initial value problem
  
  %numerical integration routine
  options=odeset('AbsTol',1e-6,'RelTol',1e-3,'InitialStep',0.005);
    
  [~,y]=ode45(@(s,y) deriv(y,k),[Break_Points(k),Break_Points(k+1)],y_init,options); 
  
  y= unique(y,'rows','stable'); %deleting double lines, 'stable': same order as before
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

  %external force- and momentdistribution
  fe=[0;0;0];
  le=[0;0;0];
    
  % Differential equations 
  v_dot=-Kse\(lie(u)*Kse*(v-[0;0;1])+R'*fe);
  u_dot=-Kbt\(lie(u)*Kbt*u+lie(v)*Kse*(v-[0;0;1])+R'*le);
  p_dot=R*v;
  R_dot=R*lie(u);
  s_dot=1;
    
  dy=[p_dot;reshape(R_dot,9,1);v_dot;u_dot;s_dot];

  end

 
end

