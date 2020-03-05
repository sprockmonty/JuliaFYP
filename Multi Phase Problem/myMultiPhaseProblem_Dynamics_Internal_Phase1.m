
function [dx,g_eq,g_neq] = myMultiPhaseProblem_Dynamics_Internal_Phase1(x,u,p,t,vdat)
% Syntax:  
%          [dx] = myProblem_Dynamics_Internal(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%    g_eq - constraint function for equality constraints
%    g_neq - constraint function for inequality constraints
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

%Constant data
mu       = vdat.mu;      %gravitational parameter
g0       = vdat.g0;      %gravity acceleration at sea level
Aref     = vdat.Aref;    %rocket reference area
rho0     = vdat.rho0;    %atmospheric density at sea level
Re       = vdat.Re;      %equatorial radius of the Earth
h0       = vdat.h0;      %density scale height
Isp      = vdat.Isp;     %specific impulse of engine in vacuum
CD       = vdat.CD;      %drag coefficient
m01      = vdat.m01;     %initial masses of first and second stage
m02      = vdat.m02;   
mfuel1   = vdat.mfuel1;  %limitive fuel masses for both stages
mfuel2   = vdat.mfuel2; 
qmax     = vdat.qmax;    %maximum dynamic pressure
accmax   = vdat.accmax;  %maximum visual acceleration
thetags  = vdat.thetags; %glide slope angle (in degrees)
rf1      = vdat.rf1;     %first stage landing site location
tmax     = vdat.tmax;    %maximum thrust
tmin     = vdat.tmin;    %minimum thrust




%Define states
%1st stage
r1 = x(:,1:3); %position vector x,y,z
v1 = x(:,4:6); %velocity vector u,v,w
m1 = x(:,7); %rocket mass 
%2nd stage
r2 = x(:,8:10); %position vector x,y,z
v2 = x(:,11:13); %velocity vector u,v,w
m2 = x(:,14); %rocket mass


%Define inputs
%1st stage
T1 = u(:,1:3); %thrust vector x force,y force,z force
%2nd stage
T2 = u(:,4:6); %thrust vector x force,y force,z force



%Define ODE right-hand side
%1st stage
h1   = rowNorm(r1) - Re;                  %altitude
rho1 = rho0*exp(-h1/h0);               %atmospheric density
D1   = -0.5*CD*Aref*rho1.*rowNorm(v1).*v1;   %drag force vector



%1st stage
rdot1 = v1;
vdot1 = -mu*r1./(rowNorm(r1)).^3 + T1./(m1+m2) + D1./(m1+m2);
mdot1 = -rowNorm(T1)/(g0*Isp);
%2nd stage
%mdot2 = -rowNorm(T2)/(g0*Isp);


dx = [rdot1,vdot1,mdot1,rdot1,vdot1,0];

%Define Path constraints
%g_eq(:,1)=g_eq1(x1,...,u1,...p,t);
%g_eq(:,2)=g_eq2(x1,...,u1,...p,t);
%...

%1st stage
ReConst1 = Re^2 - norm(r1)^2;                         %altitude constraint 
qmaxConst1 = 0.5*rho1.*rowNorm(v1).^2 - qmax;         %maximum dynamic pressure constraint
accmaxConst1 = rowNorm((T1+D1)./m1).^2 - accmax.^2;   %maximum visual acceleration constraint
thrustConst1 = (9*tmax)^2 - rowNorm(T1)^2;

%2nd stage
ReConst2 = Re^2 - rowNorm(r2)^2;                      %altitude constraint
qmaxConst2 = 0.5*rho2.*rowNorm(v2).^2 - qmax;         %maximum dynamic pressure constraint
accmaxConst2 = rowNorm((T2+D2)./m2).^2 - accmax.^2;   %maximum visual acceleration constraint

g_eq = [thrustConst1];

g_neq=[ReConst1,ReConst2,...
       qmaxConst1,qmaxConst2,...
       accmaxConst1,accmaxConst2,...
       ];

function a = rowNorm(a) 
    %takes the norm of each row in an nx3 array
    %returns array where each element is the norm of the corresponding row
    a = arrayfun(@(x,y,z) norm([x,y,z]), a(:,1),a(:,2),a(:,3));
end
%------------- END OF CODE --------------
end
