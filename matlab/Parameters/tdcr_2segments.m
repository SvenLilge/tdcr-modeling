
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga


% Geometrical and mechanical properties of the 2 sections TDCR

% Subsegment number and length
n = [10;10]; %number of spacerdisks on the 2 sections, n(1) corresponds to proximal section
n_disk=sum(n); %number of spacerdisks
L = [200e-3;200e-3]; % Lenght of the section in m
l=[L(1)/n(1)*ones(1,n(1)) L(2)/n(2)*ones(1,n(2))]; %array of segment lengths of each subsegment

% Tendon position on the disks
r_disk=10*1e-3; %distance of routing channels to backbone-center[m]
% Tendons actuating the proximal segment
r1=[r_disk 0 0 1]'; 
r2=[0 r_disk 0 1]'; 
r3=[-r_disk 0 0 1]'; 
r4=[0 -r_disk 0 1]'; 
% Tendons actuating the distal segment
r5=[r_disk*cos(pi/4) r_disk*sin(pi/4) 0 1]'; 
r6=[r_disk*cos(3*pi/4) r_disk*sin(3*pi/4) 0 1]'; 
r7=[r_disk*cos(5*pi/4) r_disk*sin(5*pi/4) 0 1]'; 
r8=[r_disk*cos(-pi/4) r_disk*sin(-pi/4) 0 1]'; 
p_tendon=[r1 r2 r3 r4 r5 r6 r7 r8]; %additional tendons can be added through additional columns

% Tendons actuating the proximal segment
r1=[0 r_disk 0 1]'; 
r2=[r_disk*cos(-pi/6) r_disk*sin(-pi/6) 0 1]';
r3=[r_disk*cos(7*pi/6) r_disk*sin(7*pi/6) 0 1]'; 
% Tendons actuating the distal segment
r4=[0 r_disk 0 1]'; 
r5=[r_disk*cos(-pi/6) r_disk*sin(-pi/6) 0 1]';
r6=[r_disk*cos(7*pi/6) r_disk*sin(7*pi/6) 0 1]';
p_tendon=[r1 r2 r3 r4 r5 r6];

% Tendon tension
% F = [8 2 8 2 0 0 0 0];
F = [8 0 0 0 0 0];

% Backbone mechanical and geometrical properties
E=54*10^9;% Youngs modulus
nu=0.3; %Poissons ratio
G=E/(2*(1+nu)); %Shear modulus
ro=1.4/2*10^-3; %outer radius of bb
ri=0; %inner radius of bb
I=1/4*pi*(ro^4-ri^4); %moment of inertia
g=0; %9.81; %acceleration due to gravity
m_bb=0.0115*l*g; %mass of backbone %weight of the backbone expressed in kg/m multiplied by length and g
m_disk=0.2*1e-3*ones(1,n_disk)*g; %array of masses of each spacerdisk

% External tip forces and moments
Ftex = [0;0;0;0]; % Force applied at the tip, expressed in global frame
Mtex = [0;0;0;0]; % Moment applied at the tip