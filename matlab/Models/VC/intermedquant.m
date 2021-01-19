
% This code implements different approaches to model the kinematics/statics
% of a two segment tendon driven continuum robot and is part of the following
% publication:

% How to model tendon-driven continuum robots and benchmark modelling performance
% Priyanka Rao, Quentin Peyron, Sven Lilge, Jessica Burgner-Kahrs
% frontiers in Robotics and AI 2021
% DOI: 10.3389/frobt.2020.630245

% Copyright (C) 2021 Continuum Robotics Laboratory, University of Toronto Mississauga

function [A,B,G,H,c,d] = intermedquant2(u,v,R,r,Kse,Kbt,tau_stell,fe,le,k)
%Calculates intermediate vector and matrix quantities depending on which
%sections are present
%----INPUT----
% k: numer of section
% u,v:angular and linear rate of change
% R: orientation
% param: parameters of manipulator
% Kse,Kbt: stiffness matrices
% refconfig: reference configuartion
% tau_stell: tendon tension
% fe,le: external force- and momentdistribution
%----OUTPUT----
% A,B,G,H: intermediate matrix quantities
% c,d: intermediate vector quantities



% tendon paths in body frame
% pib_dot=lie(u)*ri+v;
p1b_dot=[-u(3)*r(2,1);u(3)*r(1,1);-u(2)*r(1,1)+u(1)*r(2,1)]+v;
p2b_dot=[-u(3)*r(2,2);u(3)*r(1,2);-u(2)*r(1,2)+u(1)*r(2,2)]+v;
p3b_dot=[-u(3)*r(2,3);u(3)*r(1,3);-u(2)*r(1,3)+u(1)*r(2,3)]+v;

%% For all sections
% (tendons of section 2 are always present)

% Ai=-tau_i*((lie(pib_dot))^2/norm(pib_dot)^3)
A21=-(tau_stell(1,2)/norm(p1b_dot)^3)*[0 -p1b_dot(3) p1b_dot(2);p1b_dot(3) 0 -p1b_dot(1);-p1b_dot(2) p1b_dot(1) 0]*[0 -p1b_dot(3) p1b_dot(2);p1b_dot(3) 0 -p1b_dot(1);-p1b_dot(2) p1b_dot(1) 0];
A22=-(tau_stell(2,2)/norm(p2b_dot)^3)*[0 -p2b_dot(3) p2b_dot(2);p2b_dot(3) 0 -p2b_dot(1);-p2b_dot(2) p2b_dot(1) 0]*[0 -p2b_dot(3) p2b_dot(2);p2b_dot(3) 0 -p2b_dot(1);-p2b_dot(2) p2b_dot(1) 0];
A23=-(tau_stell(3,2)/norm(p3b_dot)^3)*[0 -p3b_dot(3) p3b_dot(2);p3b_dot(3) 0 -p3b_dot(1);-p3b_dot(2) p3b_dot(1) 0]*[0 -p3b_dot(3) p3b_dot(2);p3b_dot(3) 0 -p3b_dot(1);-p3b_dot(2) p3b_dot(1) 0];
A=A21+A22+A23;

% Bi=lie(ri)*Ai
B21=[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]*A21;
B22=[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]*A22;
B23=[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]*A23;
B=B21+B22+B23; 

% G= -sum Ai*lie(ri)
G =-(A21*[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]+...
   A22*[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]+...
   A23*[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]);
% H= -sum Bi*lie(ri)
H =-(B21*[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]+...
   B22*[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]+...
   B23*[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]);

% ai=Ai*lie(u)*pib_dot
a21=A21*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p1b_dot;
a22=A22*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p2b_dot;
a23=A23*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p3b_dot;
a=a21+a22+a23;

% b=lie(ri)*ai
b=[r(2,1)*a21(3);-r(1,1)*a21(3);-r(2,1)*a21(1)+r(1,1)*a21(2)]+...
    [r(2,2)*a22(3);-r(1,2)*a22(3);-r(2,2)*a22(1)+r(1,2)*a22(2)]+...
    [r(2,3)*a23(3);-r(1,3)*a23(3);-r(2,3)*a23(1)+r(1,3)*a23(2)];

%% Additions depending on current section

if k==1  %section 1
    
    % Ai=-tau_i*((lie(pib_dot))^2/norm(pib_dot)^3)
    A11=-(tau_stell(1,1)/norm(p1b_dot)^3)*[0 -p1b_dot(3) p1b_dot(2);p1b_dot(3) 0 -p1b_dot(1);-p1b_dot(2) p1b_dot(1) 0]*[0 -p1b_dot(3) p1b_dot(2);p1b_dot(3) 0 -p1b_dot(1);-p1b_dot(2) p1b_dot(1) 0];
    A12=-(tau_stell(2,1)/norm(p2b_dot)^3)*[0 -p2b_dot(3) p2b_dot(2);p2b_dot(3) 0 -p2b_dot(1);-p2b_dot(2) p2b_dot(1) 0]*[0 -p2b_dot(3) p2b_dot(2);p2b_dot(3) 0 -p2b_dot(1);-p2b_dot(2) p2b_dot(1) 0];
    A13=-(tau_stell(3,1)/norm(p3b_dot)^3)*[0 -p3b_dot(3) p3b_dot(2);p3b_dot(3) 0 -p3b_dot(1);-p3b_dot(2) p3b_dot(1) 0]*[0 -p3b_dot(3) p3b_dot(2);p3b_dot(3) 0 -p3b_dot(1);-p3b_dot(2) p3b_dot(1) 0];
    A=A+A11+A12+A13;
    
    % Bi=lie(ri)*Ai
    B11=[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]*A11;
    B12=[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]*A12;
    B13=[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]*A13;
    B=B+B11+B12+B13;
    
    % G= -sum Ai*lie(ri)
    G =G-(A11*[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]+...
        A12*[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]+...
        A13*[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]);
    
    % H= -sum Bi*lie(ri)
    H =H-(B11*[0 0 r(2,1);0 0 -r(1,1);-r(2,1) r(1,1) 0]+...
        B12*[0 0 r(2,2);0 0 -r(1,2);-r(2,2) r(1,2) 0]+...
        B13*[0 0 r(2,3);0 0 -r(1,3);-r(2,3) r(1,3) 0]);
    
    % ai=Ai*lie(u)*pib_dot
    a11=A11*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p1b_dot;
    a12=A12*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p2b_dot;
    a13=A13*[0 -u(3) u(2);u(3) 0 -u(1); -u(2) u(1) 0]*p3b_dot;
    a=a+a11+a12+a13;
   
    % b=lie(ri)*ai
    b=b+[r(2,1)*a11(3);-r(1,1)*a11(3);-r(2,1)*a11(1)+r(1,1)*a11(2)]+...
        [r(2,2)*a12(3);-r(1,2)*a12(3);-r(2,2)*a12(1)+r(1,2)*a12(2)]+...
        [r(2,3)*a13(3);-r(1,3)*a13(3);-r(2,3)*a13(1)+r(1,3)*a13(2)];
    
end

% v = [0;0;1];
c=-lie(u)*Kbt*u-lie(v)*Kse*(v-[0;0;1])-R'*le-b;
d=-lie(u)*Kse*(v-[0;0;1])-R'*fe-a;


end