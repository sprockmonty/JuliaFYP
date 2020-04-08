function [problem,guess] = myMultiPhaseProblem_Phase2(problem_mp, guess_mp)
%myProblem - Template file for optimal control problem definition
%
%Syntax:  [problem,guess] = myProblem
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% MAT-files required: none
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

%Final varibles
r1f        = [5605.2e3,0,3043.4e3];
v1f        = [0, 0, 0];
af2        = 6593.145e3;
ef2        = 0.0076;
if2        = 28.5;
rightAscf2 = 269.8;
omegaf2    = 130.5;
m10 = 431.6e3;
m20 = 107.5e3;
mfuel1 = problem_mp.data.mfuel1;
mfuel2 = problem_mp.data.mfuel2;


% Plant model name, provide in the format of function handle
InternalDynamics=@myMultiPhaseProblem_Dynamics_Internal_Phase2;
SimDynamics=[];

% Analytic derivative files (optional), provide in the format of function handle
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_myMultiPhaseProblem_Phase2;

% Initial time. t0<tf
problem.time.t0_idx=2;
problem.time.t0_min=problem_mp.time.t_min(problem.time.t0_idx);
problem.time.t0_max=problem_mp.time.t_max(problem.time.t0_idx); %separation could occur at any point between start and end time
guess.t0=guess_mp.time(problem.time.t0_idx);

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_idx=3;
problem.time.tf_min=problem_mp.time.t_min(problem.time.tf_idx);     
problem.time.tf_max=problem_mp.time.t_max(problem.time.tf_idx); 
guess.tf=guess_mp.time(problem.time.tf_idx);

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[-inf,-inf,-inf,-inf,-inf,-inf,m10 - mfuel1,-inf,-inf,-inf,-inf,-inf,-inf,m20 - mfuel2]; 
problem.states.x0u=[inf,inf,inf,inf,inf,inf,m10,inf,inf,inf,inf,inf,inf,m20]; 


% State bounds. xl=< x <=xu
problem.states.xl=[-inf,-inf,-inf,-inf,-inf,-inf,m10 - mfuel1,-inf,-inf,-inf,-inf,-inf,-inf,m20 - mfuel2];
problem.states.xu=[inf,inf,inf,inf,inf,inf,m10,inf,inf,inf,inf,inf,inf,m20];

% State rate bounds. xrl=< x <=xru
% problem.states.xrl=[]; 
% problem.states.xru=[]; 

% State error bounds
problem.states.xErrorTol_local=ones(1,14)*1e-2; 
problem.states.xErrorTol_integral=ones(1,14)*1e-2; 

% State constraint error bounds
problem.states.xConstraintTol=ones(1,14)*1e-2;
% problem.states.xrConstraintTol=[eps_x1dot_bounds ... eps_xndot_bounds];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[r1f,v1f,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf]; 
problem.states.xfu=[r1f,v1f,inf,inf,inf,inf,inf,inf,inf,inf];

% Guess the state trajectories with [x0 ... xf]
% guess.time=[t0 ... tf];
guess.states(:,1)=[r1f(1) r1f(1)];
guess.states(:,2)=[r1f(2) r1f(2)];
guess.states(:,3)=[r1f(3) r1f(3)];
guess.states(:,4)=[v1f(1) v1f(1)];
guess.states(:,5)=[v1f(2) v1f(2)];
guess.states(:,6)=[v1f(3) v1f(3)];
guess.states(:,7)=[m10-mfuel1 m10-mfuel1];
guess.states(:,8)=[r1f(1) r1f(1)];
guess.states(:,9)=[r1f(2) r1f(2)];
guess.states(:,10)=[r1f(3) r1f(3)];
guess.states(:,11)=[v1f(1) v1f(1)];
guess.states(:,12)=[v1f(2) v1f(2)];
guess.states(:,13)=[v1f(3) v1f(3)];
guess.states(:,14)=[m20-mfuel2 m20-mfuel2];


% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[0,0,0,0,0,0];
problem.inputs.uu=[inf,inf,inf,inf,inf,inf];

% Bounds on the first control action
problem.inputs.u0l=[-inf,-inf,-inf,-inf,-inf,-inf];
problem.inputs.u0u=[inf,inf,inf,inf,inf,inf];

% Input rate bounds
% problem.inputs.url=[u1dot_lowerbound ... umdot_lowerbound]; 
% problem.inputs.uru=[u1dot_upperbound ... umdot_upperbound]; 

% Input constraint error bounds
problem.inputs.uConstraintTol=ones(1,6)*1e-2;
% problem.inputs.urConstraintTol=[eps_u1dot_bounds ... eps_umdot_bounds];

% Guess the input sequences with [u0 ... uf]
guess.inputs(:,1)=[360e3 360e3];
guess.inputs(:,2)=[360e3 360e3];
guess.inputs(:,3)=[360e3 360e3];
guess.inputs(:,4)=[934e3 934e3];
guess.inputs(:,5)=[934e3 934e3];
guess.inputs(:,6)=[934e3 934e3];

% Path constraint function 
problem.constraints.ng_eq=1; % number of quality constraints in format of g(x,u,p,t) == 0
problem.constraints.gTol_eq=[1e-2]; % equality cosntraint error bounds

problem.constraints.gl=[-inf,-inf,-inf,-inf,-inf,-inf,-inf,problem_mp.data.tmin^2]; % Lower ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gu=[0,0,0,0,0,0,0,problem_mp.data.tmax^2]; % Upper ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gTol_neq=ones(1,6)*1e-2; % inequality constraint error bounds

% OPTIONAL: define the time duration each constraint will be active, for
% example (for ECH enabled in setings)
%problem.constraints.gActiveTime{1}=[guess.tf/2 guess.tf];
%problem.constraints.gActiveTime{2}=[];
%...
%problem.constraints.gActiveTime{5}=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[af2,ef2,if2,rightAscf2,omegaf2];
problem.constraints.bu=[af2,ef2,if2,rightAscf2,omegaf2];
problem.constraints.bTol=[1e-2,1e-2,1e-2,1e-2,1e-2]; 

% store the necessary problem parameters used in the functions
problem.data=problem_mp.data;
% optional setting for automatic regularization
problem.data.penalty.values=[];
problem.data.penalty.i=1; %starting weight

% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.data.functions_unscaled=problem.functions_unscaled;
problem.data.ng_eq=problem.constraints.ng_eq;
problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];


%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,vdat)

% L_unscaled - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;

%------------- BEGIN CODE --------------

%Define states and setpoints
%x1 = x(:,1); 
% %...
%xn=x(:,n); 
%
% %Define inputs
%u1 = u(:,1);
% % ...
%um = u(:,m);

stageCost = 0*t;

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 

% E_unscaled - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

boundaryCost=-xf(14);

%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
% Leave it here
varargin=varargin{1};
%------------- BEGIN CODE --------------
[a,e,i,rightAsc,omega] = earth2orbital(xf(8:10),xf(11:13));
bc(1,:)=a;
bc(2,:)=e;
bc(3,:)=i;
bc(4,:)=rightAsc;
bc(5,:)=omega;

if any(any(isnan(bc)))
    pause;
end
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if strcmp(options.transcription,'hpLGR') && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

