
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

clear all;
close all;

addpath('Parameters/');
addpath('Models/CCsub/');
addpath('Models/VC/');
addpath('Models/PRBM/');
addpath('Models/CC/');

%% Load TDCR parameters
tdcr_2segments;

% Tendon forces
F = [2 0 0 0 1 0];

% External forces and moments at the tip
Ftex = [0;0;0;0];
Mtex = [0;0;0;0];

%% Solving the CCsub model
var0 = 0.01*ones(1,3*n_disk);
var_ccsub = CCsub_solver(n,l,F,p_tendon,m_disk,m_bb,E,I,G,Ftex,Mtex,var0);

%% Plotting
figure(1);
plot_tdcr_ccsub(l,p_tendon(1:3,:),n_disk,var_ccsub,r_disk);
view(90,0);

% Construct the robot shape with CCsub output, using 30 intermediate
% points along each subsegment
s_disk = [0 L(1)/n(1):L(1)/(n(1)):L(1) L(1)+L(2)/n(2):L(2)/n(2):sum(L)];
[stest,pos_ccsub] = construct_tdcr_ccsub(l,p_tendon(1:3,:),n_disk,var_ccsub);

%% Solving the VC models

% Set up TDCR parameters
param=setup_param_vc(0,0,0,0,sum(L),ro,ri,m_bb(1),E,nu,0,r_disk,m_disk(1),n_disk,L,r1(1:3),r2(1:3),r3(1:3));
Fc = [F(1:3);F(4:6)]';

% Shooting method for model solution

%initial guess for linear and angular rate of change
v_init_guess=[0;0;1]; 
u_init_guess=[0;0;0]; 
init_guess=[v_init_guess;u_init_guess];

% Solving the VC model
[y,vu_f] = VC_solver(init_guess,Fc,Ftex(1:3),L,param);
sc_vc= y(:,19); 
g_vc = [y(:,4:6) zeros(size(y,1),1) y(:,7:9) zeros(size(y,1),1) y(:,10:12) zeros(size(y,1),1) y(:,1:3) ones(size(y,1),1)];

% Solving the VCref model
init_guess = repmat(init_guess,sum(n),1);
[y,vu_f] = VCref_solver(init_guess,F,Ftex(1:3),Mtex(1:3),L,param,n);
sc_vcref= y(:,19); 
g_vcref = [y(:,4:6) zeros(size(y,1),1) y(:,7:9) zeros(size(y,1),1) y(:,10:12) zeros(size(y,1),1) y(:,1:3) ones(size(y,1),1)];


%% Plotting
figure(2);
plot_tdcr_vc(g_vc,sc_vc,L)
view(90,0);

figure(3);
plot_tdcr_vc(g_vcref,sc_vcref,L)
view(90,0);

%% Solving the PRBM model

% Number of rigid bodies
nrb = 4;
% Length of rigid bodies, optimized with particle swarm algorithm in Chen
% 2011
gamma= [0.125 0.35 0.388 0.136]/sum([0.136 0.388 0.35 0.125]);

%phi0 = pi/2;
phi0 = 0;

% Initialization
var0 = repmat([phi0 0*ones(1,nrb)],1,n_disk);

[var_prbm,exitflag,res] = PRBM_solver(n,nrb,gamma,l,F,p_tendon,m_disk,m_bb,E,I,G,Ftex,Mtex,var0);

%% Plotting
figure(4);
plot_tdcr_prbm(nrb,gamma,l,p_tendon,n_disk,var_prbm,r_disk);

% Compute the position of each disk and each node of the PRB model
[stest,pos_prbm] = construct_tdcr_prbm(nrb,gamma,l,n_disk,var_prbm);


%% Solving the constant curvature model (geometric)

% Compute tendon length from the robot shape obtained with Cosserat rod
% model
[delta_l] = calctendonlength(g_vc,sc_vc,param,L);
l_d = [-delta_l(1:3)/1000+L(1);-(delta_l(4:6)-delta_l(1:3))/1000+L(2)];

var_cc = CC_solver(l_d,n,r_disk);

% Construct the robot shape with CC output, using 50 points along each
% segment
[T1_cc,T2_cc] = construct_tdcr_cc(var_cc);

%% Plotting
figure(5);
plot_tdcr_cc(T1_cc,T2_cc,L);



