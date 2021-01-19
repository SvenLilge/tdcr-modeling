
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [p_f] = plot_tdcr_prbm(nrb,gamma,l,p_tendon,n_disk,var,r_disk)

%% visualisation of the results
n_tendon = size(p_tendon,2);
p_last=[0;0;0];
[T_i,Trb]=trans_mat_prbm(var,nrb,gamma,l,n_disk,0);
p_plot = p_last;
for j=1:size(Trb,3)
    p_plot = [p_plot Trb(1:3,4,j)];
end

p_f = p_last;
for ss_i=1:n_disk
    
    % Tendon initial and final position
    T_i=trans_mat_prbm(var,nrb,gamma,l,ss_i-1,0); %returns transformation matrix from i to i-1
    p_t1=T_i*p_tendon; % coordinates of tends on disk ss_i-1
    T_i1=trans_mat_prbm(var,nrb,gamma,l,ss_i,0);
    p_t2=T_i1*p_tendon;

    
    plot3(p_plot(1,:),p_plot(2,:),p_plot(3,:),'k','LineWidth',3) %plotting the backbone
    hold on;
        color = [1 1 1]*0.9;
        squaresize = 0.07;
        lmax=sum(l);
        axis([-lmax lmax -lmax lmax -lmax*2/3 4/3*lmax])
        fill3([1 1 -1 -1]*squaresize,[-1 1 1 -1]*squaresize,[0 0 0 0],color);
        xlabel('x (m)')
        ylabel('y (m)')
        zlabel('z (m)')
        grid on
        view([1 1 1])
        daspect([1 1 1])  
    p_last=p_plot(:,nrb*ss_i+1);
    p_f = [p_f p_plot];
    
    %plotting disks
    t_disk=0:0.01:2*pi; %parameter for plotting disk i
    p_disk=repmat(p_last,1,length(t_disk))+(r_disk+0.002)*T_i1(1:3,1:3)*[cos(t_disk);sin(t_disk);zeros(size(t_disk))];
    patch(p_disk(1,:)',p_disk(2,:)',p_disk(3,:)',[0 0.8 0.8]);

    %plotting tendons
    plot3([p_t1(1,:);p_t2(1,:)],[p_t1(2,:);p_t2(2,:)],[p_t1(3,:);p_t2(3,:)],'r','LineWidth',1)
end
end