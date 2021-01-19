
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [p_plot,p] = plot_tdcr_ccsub(l,p_tendon,n_disk,var,r_disk)

%% visualisation of the results
n_tendon = size(p_tendon,2);
p_last=[0;0;0];
p = [0;0;0];
for ss_i=1:n_disk
    beta=var(3*ss_i-2);
    gamma=var(3*ss_i-1);
    k=sqrt(beta^2+gamma^2);
    phi=atan2(gamma,beta);
    r=1/k;
    s=0.0:0.001:l(ss_i);
    p_i=r*[(1-cos(s/r))*cos(phi);(1-cos(s/r))*sin(phi);sin(s/r)]; %coordinates of points along subsegment i
    R_i=rotation_mat_ccsub(var,l,ss_i-1,0);
    p_plot=R_i*p_i+repmat(p_last,1,length(s)); % coordinates of end of ss_i
    p_t1=R_i*p_tendon+repmat(p_last,1,n_tendon); % coordinates of tends on disk ss_i-1
    
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
    p_last=p_plot(:,end);
    p = [p p_plot(:,end)];
    
    %plotting disks
    R_i1=rotation_mat_ccsub(var,l,ss_i,0);
    t_disk=0:0.01:2*pi; %parameter for plotting disk i
    p_disk=repmat(p_last,1,length(t_disk))+(r_disk+0.01)*R_i1*[cos(t_disk);sin(t_disk);zeros(size(t_disk))];
    patch(p_disk(1,:)',p_disk(2,:)',p_disk(3,:)',[0 0.8 0.8]);
    p_t2=R_i1*p_tendon+repmat(p_last,1,n_tendon); % coordinates of tendons on disk ss_i
    
    %plotting tendons
    plot3([p_t1(1,:);p_t2(1,:)],[p_t1(2,:);p_t2(2,:)],[p_t1(3,:);p_t2(3,:)],'r','LineWidth',1)
end
end