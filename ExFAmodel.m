% Here we run experiments for the SIREV paper. This is the FA model run on
% a degree regular network (z=4), with the purpose being to compare the
% AME, PA and MF frameworks. 

% MAKE SURE TO UNCOMMENT LINES 20-27 IN THE F_rates.m FILE

% Get the degree distribution Pk aswell as Kmin, Kmaz, average degree etc
DegreeDistribution = 'zRegular'; % 'PRG' or 'truncSFN' or 'zRegular' or 'custom'
z = 4;
DistParams = [z];   % Mean Degree

% Parameters
f = 2;
T = 0.4;
mag0 = 1/(1+exp(-1/T));
rho0 = [1-mag0,mag0,0,0]; 

% Simulation inputs
n=4;
endtime = power(10,6);

% ---------------
% EXPERIMENTS
% ---------------

% Mean field
[TMF xMF] = multi_state_solver(n,DegreeDistribution,z,[T f], rho0, endtime, 'MF');
PhiMF = xMF(:,1) + xMF(:,2);

% AME
[TAME xAME] = multi_state_solver(n,DegreeDistribution,z,[T f], rho0, endtime, 'AME');
PhiAME = xAME(:,1) + xAME(:,2);

% PA
% Need to change ICs to have no zeros (as messes with PA)
eps = 0.0001;
rho0 = [1-mag0-eps,mag0-eps,eps,eps]; 
[TPA xPA] = multi_state_solver(n,DegreeDistribution,z,[T f], rho0, endtime, 'PA');
PhiPA = xPA(:,1) + xPA(:,2);

% MC simulations
load('simulation_outputs/phi.txt')
phi = phi(2:length(phi));
tspan = logspace(0, log10(length(phi)), 25);
entries = unique(round(tspan));
PhiMC = phi(entries);
PhiMC = [PhiMC; PhiMC(length(PhiMC))*ones(5,1)];
entries = [entries logspace(log10(length(phi)),log10(length(phi)*10), 5)];

% Mean Field plot
figure
semilogx(TAME,PhiAME,'r-',TPA,PhiPA,'b-.',TMF,PhiMF,'g--',entries, PhiMC,'kx','LineWidth',1)
xlabel('$t$','Interpreter','Latex');
ylabel('$\phi(t)\;\;\;$  ','Interpreter','Latex')
ylim([0 1])
l = legend('AME','PA','MF','Simulation');
set(l,'Interpreter','latex','Location','southwest')
legend boxoff
set(get(gca,'ylabel'),'rotation',0)


