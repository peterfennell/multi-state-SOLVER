% -------------------------------------------------------------------------
% This is the Cooperative SIS model as specified in Section 4.1 of [1].
% Inputs to the model are the base infectivity rate beta and the
% accentutation parameter lambda

% PLEASE MAKE SURE THE CORRECT RATE FUNCTION IS SPECIFIED IN THE FILE F_rates.m
% BY COMMENTING OUT LINES 43-48

% [1] Fennell, Gleeson (2017) "MULTISTATE DYNAMICAL PROCESSES ON NETWORKS:
% ANALYSIS THROUGH DEGREE-BASED APPROXIMATION FRAMEWORKS" 
% -------------------------------------------------------------------------

% Add the path of the functions to the search path
addpath ../multi-state-SOLVER

% Get the degree distribution Pk aswell as Kmin, Kmaz, average degree etc
DegreeDistribution = 'zRegular'; % 'PRG' or 'truncSFN' or 'zRegular' or 'custom'
z = 4;
DistParams = [z];   % Mean Degree

% % Co-SIS parameters
beta = 0.9/z;
lambda = 5;
DynamicsParams = [beta lambda];

% unstable equilibrium
iminus = (lambda-2)./(2*(lambda-1)) - sqrt((z*lambda*beta)^2-4*(lambda-1))./(2*z*beta*(lambda-1));

% Simulation inputs
n=4;
endtime = 100;

% ---------------
% EXPERIMENTS
% ---------------

% 1: At criticality
s0 = 1-iminus;
b0 = (1-s0)*(1-z*s0*beta)/(1+z*s0*beta);
x10 = (iminus-b0)/2;
x20 = (iminus-b0)/2;
rho0 = [s0 x10 x20 b0];

% ODE SOLVER
[TMFcrit xMFcrit] = multi_state_solver(n,DegreeDistribution,z,[beta lambda], rho0, endtime, 'MF');
IMFcrit = 1-xMFcrit(:,1);

% 2: Above criticality
gamma = 1.01;
s0 = 1-(gamma*iminus);
b0 = (1-s0)*(1-z*s0*beta)/(1+z*s0*beta);
x10 = (gamma*iminus-b0)/2;
x20 = (gamma*iminus-b0)/2;
rho0 = [s0 x10 x20 b0];

% ODE SOLVER
[TMFsuper xMFsuper] = multi_state_solver(n,DegreeDistribution,z,[beta lambda], rho0, endtime, 'MF');
IMFsuper = 1-xMFsuper(:,1);

% 3: Below criticality
gamma = 0.99;
s0 = 1-(gamma*iminus);
b0 = (1-s0)*(1-z*s0*beta)/(1+z*s0*beta)
x10 = (gamma*iminus-b0)/2;
x20 = (gamma*iminus-b0)/2;
rho0 = [s0 x10 x20 b0];

% ODE SOLVER
[TMFsub xMFsub] = multi_state_solver(n,DegreeDistribution,z,[beta lambda], rho0, endtime, 'MF');
IMFsub = 1-xMFsub(:,1);


% Mean Field plot
figure
plot(TMFsub,IMFsub,'k-.', TMFcrit,IMFcrit,'k-',TMFsuper,IMFsuper,'k--','LineWidth',1)
xlabel('$t$','Interpreter','Latex');
ylabel('$i(t)\;\;\;$  ','Interpreter','Latex')
ylim([0 1])
l = legend('$i(0)<\bar{i}_i$','$i(0)=\bar{i}_i$','$i(0)>\bar{i}_i$')
set(l,'Interpreter','latex','Location','northwest')
set(get(gca,'ylabel'),'rotation',0)

