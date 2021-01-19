
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [var,exitflag] = CCsub_solver(n,l,F,p_tendon,m_disk,m_bb,E,I,G,Ftex,Mtex,var0)
%%returns the solved values for [beta, gamma, epsi] for each subsegment i
n_disk = sum(n);

options1 = optimset('Display','iter','TolFun',1e-6); 

tic
[var,res,exitflag,output]=fsolve(@optim_f,var0,options1);
toc

%% solver

    function [res] = optim_f(var)
        res=zeros(n_disk*3,1);
        for ss_i=n_disk:-1:1 %iterating over each subsegment
            
            % Kinematics
            T_i=trans_mat_ccsub(var,l,ss_i,ss_i-1); %returns transformation matrix rotation from i to i-1
            %initialising arc parameters as defined in [2]
            beta=var(3*ss_i-2);
            gamma=var(3*ss_i-1);
            k=sqrt(beta^2+gamma^2);
            phi=atan2(gamma,beta);
            theta=l(ss_i)*k;
            epsi=var(3*ss_i);
            
            p_ti=T_i*p_tendon; %position of tendon k at diski wrt i-1 frame
            p_i=T_i*[0;0;0;1]; %position of diski wrt i-1 frame
            norm_ct1=sqrt(sum((-p_ti(1:3,:)+p_tendon(1:3,:)).^2));
            
            % Tendon tension on the different disks
            % Tension of the first 4 tendons apply on the proximal segment
            % only
            nt = length(F)/length(n); % Number of tendons per section
            Fdisk = [repmat(F,n(1),1);[zeros(n(2),nt) repmat(F(nt+1:end),n(2),1)]];  
            
            % Direction orthogonal to the disk
            zi= T_i(:,3);
            
            if ss_i<n_disk
                T_i1=trans_mat_ccsub(var,l,ss_i+1,ss_i-1); %rotation from i+1 to i-1, to express all the forces in the same frame (of disk i-1)
                p_ti1=T_i1*p_tendon;%position of tendon k at diski+1
                norm_ct2=sqrt(sum((-p_ti(1:3,:)+p_ti1(1:3,:)).^2));
                
                F_rel=((-p_ti+p_tendon)./repmat(norm_ct1,4,1)).*repmat(Fdisk(ss_i,:),4,1)+((-p_ti+p_ti1)./repmat(norm_ct2,4,1)).*repmat(Fdisk(ss_i+1,:),4,1);%.*repmat(F,4,1); %net force due to tendons on disk i
                if ss_i==n(1)
                    % Tip of segment 1
                    % Consider the full force for tendon 1 to 3, remove
                    % component orthogonal to the disk for tendon 4 to 6
                    F_rel = [F_rel(:,1:3) F_rel(:,4:6)-repmat((zi'*F_rel(:,4:6)),[4,1]).*repmat(zi,[1,3])];
                else
                    % Remove component orthogonal to the disk for tendon 1
                    % to 6
                    F_rel = F_rel - repmat(zi'*F_rel,[4,1]).*repmat(zi,[1,6]);
                end
            elseif ss_i==n_disk
                % Tip of segment 2
                % Consider full force of tendon 4 to 6 (Fdisk for tendon 1
                % to 3 equals 0 here)
                F_rel=((-p_ti+p_tendon)./repmat(norm_ct1,4,1)).*repmat(Fdisk(ss_i,:),4,1); %net force due to tendon 1 on disk i
                F_prev=[0;0;0;0];M_prev=[0;0;0]; %no recursive force from following subsegments
            end
            
            M_rel = cross_product(p_ti(1:3,:),F_rel(1:3,:));% moment from the tendon forces
        
            % External forces and moments
            if ss_i==n_disk
                Rt = trans_mat_ccsub(var,l,ss_i,0);
                Fex = Rt\Ftex;
                Mex = Rt(1:3,1:3)\Mtex(1:3)+cross_product(p_i(1:3),Fex(1:3));
            else
                Fex = zeros(4,1); % Apart from gravity, no external froces applied along the robot backbone
                Mex = zeros(3,1);
            end
            
            % Gravity forces and moment            
            Rg=trans_mat_ccsub(var,l,ss_i-1,0); %returns transformation matrix rotation from i-1 to 0 (ground frame)
            Fg_disk=(m_disk(ss_i))*(Rg\[0;0;-1;0]); %gravity forces due to disk i
%             Fg_disk = zeros(4,1);
            
            g_cog=sin(theta/2)/(theta/2);
            p_g=(1/k)*[(1-g_cog*cos(theta/2))*cos(phi);(1-g_cog*cos(theta/2))*sin(phi);g_cog*sin(theta/2)]; %center of gravity of backbone between two disks
            Fg_bb=(m_bb(ss_i))*(Rg\[0;0;-1;0]); %gravity forces due to part of backbone in subsegment i
%             p_g = [0;0;0;1];
%             Fg_bb = zeros(4,1);
            
            % Forces and moments of subsequent disks
            F_prev=T_i*F_prev; % transforming recursie force from subsequent disks expressed wrt i
            M_prev=T_i*[M_prev;0];  % transforming recursie moment from subsequent disks expressed wrt i
         
            % Net forces and moment due to tendon tension and external
            % forces
            F_net=sum(F_rel,2)+F_prev+Fg_disk+Fg_bb+Fex;
            M_net=sum(M_rel,2)+M_prev(1:3)+cross_product(p_i(1:3),F_prev(1:3))+cross_product(p_i(1:3),Fg_disk(1:3))+cross_product(p_g(1:3),Fg_bb(1:3))+Mex;
            
            % Moments due to backbone elasticity
            M_bend=[cos(phi) -sin(phi) 0;
                sin(phi) cos(phi) 0;
                0 0 1]*[0 k*E*I 0]'; %eqn 24 from [2]
            M_tor=T_i*[0 0 2*I*G*epsi/l(ss_i) 0]'; %eqn 25 from [2]
            
            % Static equilibrium
            res(3*ss_i-2:3*ss_i)=M_bend+M_tor(1:3)-M_net;
            F_prev=F_net;
            M_prev=M_net;

        end
    end
end