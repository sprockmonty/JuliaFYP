function [problem,guess,phaseoptions] = myMultiPhaseProblem
%myMultiPhaseProblem - Template file for optimal control problem definition using multi-phase formulation
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%    phaseoptions - options for each phases
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


% Initial and final time for different phases. Let t_min(end)=t_max(end) if tf is fixed.
problem.mp.time.t_min=[0,0,0];     
problem.mp.time.t_max=[0 700 700];  
guess.mp.time=[0 150 600];

% Parameters bounds. pl=< p <=pu
problem.mp.parameters.pl=[];
problem.mp.parameters.pu=[];
guess.mp.parameters=[];

% Bounds for linkage boundary constraints bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
problem.mp.constraints.bll.linear=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];
problem.mp.constraints.blu.linear=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];
problem.mp.constraints.blTol.linear=ones(1,14)*1e-2; 

 problem.mp.constraints.bll.nonlinear=[];
 problem.mp.constraints.blu.nonlinear=[];
 problem.mp.constraints.blTol.nonlinear=[]; 

% Get function handles
problem.mp.linkfunctions=@bclink;

% Store the necessary problem parameters used in the functions
auxdata.mu = 3.986012e14;      %gravitational parameter
auxdata.g0 = 9.801;      %gravity acceleration at sea level
auxdata.Aref = 4*pi;    %rocket reference area
auxdata.rho0 = 1.225;    %atmospheric density at sea level
auxdata.Re = 6378.145e3;      %equatorial radius of the Earth
auxdata.h0 = 7200;      %density scale height
auxdata.Isp = 340;     %specific impulse of engine in vacuum
auxdata.CD = 0.5;      %drag coefficient
auxdata.m01 = 431.6;     %initial masses of first and second stage
auxdata.m02 = 107.5;   
auxdata.mfuel1 = 409.5;  %limitive fuel masses for both stages
auxdata.mfuel2 = 103.5; 
auxdata.qmax = 80e3;    %maximum dynamic pressure
auxdata.accmax = 10 * auxdata.g0;  %maximum visual acceleration
auxdata.thetags = 80; %glide slope angle
auxdata.rf1 = [5605.2e3,0,3043.4e3];     %first stage landing site location
auxdata.tmax = 934e3;    %maximum thrust
auxdata.tmin = 360e3;    %minimum thrust
problem.mp.data = auxdata;

% Define different phases of OCP
[problem.phases{1},guess.phases{1}] = myMultiPhaseProblem_Phase1(problem.mp, guess.mp);
[problem.phases{2},guess.phases{2}] = myMultiPhaseProblem_Phase2(problem.mp, guess.mp);
...

% Each phase could use different discretization method
phaseoptions{1}=settings_myMultiPhaseProblem_Phase1(3); % h method for phase 1
phaseoptions{2}=settings_myMultiPhaseProblem_Phase2(3); % h method for phase 1



...
%------------- END OF CODE --------------


function [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)

% bclink - Returns the evaluation of the linkage boundary constraints: bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
%
% Syntax:  [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    vdat- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    blc_linear - column vector containing the evaluation of the linear linkage boundary constraint functions
%    blc_nonlinear - column vector containing the evaluation of the nonlinear linkage boundary constraint functions
%
%------------- BEGIN CODE --------------

% Variable of different phase could be called using syntex of, for example, xf{n_phase}(n): the value of nth state at tf of phase number 'n_phase'

% linear linkage constraints 
blc_linear(1,:)=xf{1}(1)-xf{2}(1);
blc_linear(2,:)=xf{1}(2)-xf{2}(2);
blc_linear(3,:)=xf{1}(3)-xf{2}(3);
blc_linear(4,:)=xf{1}(4)-xf{2}(4);
blc_linear(5,:)=xf{1}(5)-xf{2}(5);
blc_linear(6,:)=xf{1}(6)-xf{2}(6);
blc_linear(7,:)=xf{1}(7)-xf{2}(7);
blc_linear(8,:)=xf{1}(8)-xf{2}(8);
blc_linear(9,:)=xf{1}(9)-xf{2}(9);
blc_linear(10,:)=xf{1}(10)-xf{2}(10);
blc_linear(11,:)=xf{1}(11)-xf{2}(11);
blc_linear(12,:)=xf{1}(12)-xf{2}(12);
blc_linear(13,:)=xf{1}(13)-xf{2}(13);
blc_linear(14,:)=xf{1}(14)-xf{2}(14);

blc_nonlinear=[];
% nonlinear linkage constraints 
% blc_nonlinear(1,:)=blc_nonlinear_1(x0,xf,u0,uf,p,t0,tf,vdat);
% blc_nonlinear(2,:)=blc_nonlinear_2(x0,xf,u0,uf,p,t0,tf,vdat);
% ...
%------------- END OF CODE --------------
