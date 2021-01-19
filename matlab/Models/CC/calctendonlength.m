
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga


function [delta_l] = calctendonlength(g,s,param,l)
%Calculates tendon lengths (spacer disk location neglected)

%tendon location
r1=[param.tendon.r1; 1];
r2=[param.tendon.r2; 1];
r3=[param.tendon.r3; 1];

%total cable length
l1=l(1);
l2=sum(l(1:2));

%indices
si=[1:1:length(s)];
end1=find(abs(s-l1)<eps);

%Initialize lengths
l11=0;
l12=0;
l13=0;
l21=0;
l22=0;
l23=0;

for i=2:length(si)
  gi=reshape(g(si(i),:),4,4);
  gi_min=reshape(g(si(i-1),:),4,4);
  temp1=norm(gi*r1-gi_min*r1);
  temp2=norm(gi*r2-gi_min*r2);
  temp3=norm(gi*r3-gi_min*r3);
  l21=l21+temp1;
  l22=l22+temp2;
  l23=l23+temp3;
  if i<=end1(1)
    l11=l11+temp1;
    l12=l12+temp2;
    l13=l13+temp3;
  end
end

%tendon displacement: positive - shortening, negative -elongation
delta_l=[l1-l11,l1-l12,l1-l13,l2-l21,l2-l22,l2-l23]*10^3; %in mm
ltendon=struct('delta_l11',delta_l(1),'delta_l12',delta_l(2),'delta_l13',delta_l(3),...
               'delta_l21',delta_l(4),'delta_l22',delta_l(5),'delta_l23',delta_l(6));
end

