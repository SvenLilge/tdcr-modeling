
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [param] = setup_param_vc(L1,ro1,ri1,m1,L2,ro2,ri2,m2,E,nu,roh,r_disk,m_disk,n_disk,l,r1,r2,r3)
%Sorts parameters in a struct
% E=60*10^9;%60*10^9; % [Pa] 54?
% nu=0.3;
% roh=6450; %[kg/m^3]

G=E/(2*(1+nu));

%second moments of area I=Ixx=Iyy
I1=1/4*pi*(ro1^4-ri1^4); 
I2=1/4*pi*(ro2^4-ri2^4); 

%areas of cross section
A1=pi*(ro1^2-ri1^2);
A2=pi*(ro2^2-ri2^2);

%stiffness matrices
Kbt1=diag([E*I1,E*I1,G*2*I1]); 
Kbt2=diag([E*I2,E*I2,G*2*I2]);
Kse1=diag([G*A1,G*A1,E*A1]);
Kse2=diag([G*A2,G*A2,E*A2]);
% Kse1 = zeros(3);
% Kse2 = zeros(3);


tube1=struct('L',L1,'ro',ro1,'ri',ri1,'A',A1,'Kbt',Kbt1,'Kse',Kse1,'q',m1*9.81);
tube2=struct('L',L2,'ro',ro2,'ri',ri2,'A',A2,'Kbt',Kbt2,'Kse',Kse2,'q',m2*9.81);

%position of tendon ri=[xi yi 0]'
% r1=[0 r_disk 0]';
% r2=[cos(30*pi/180)*r_disk -sin(30*pi/180)*r_disk 0]';
% r3=[-cos(30*pi/180)*r_disk -sin(30*pi/180)*r_disk 0]';

tendon=struct('r1',r1(1:3),'r2',r2(1:3),'r3',r3(1:3));
    
param = struct('tube1',tube1,'tube2',tube2,'tendon',tendon,'roh',roh,'q_disks',(n_disk*m_disk*9.81)*ones(size(l))./l);
end

