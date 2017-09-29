% Running the transient dynamics for the AME, PA and MF frameworks for two different values of the coopertative SIS transmission parameter beta (one below the critical point, one above)

% MAKE SURE THAT LINES 12-17 OF F_rates.m ARE UNCOMMENTED

% Get the degree distribution Pk aswell as Kmin, Kmaz, average degree etc
DegreeDistribution = 'zRegular'; % 'PRG' or 'truncSFN' or 'zRegular' or 'custom'
z = 8;
DistParams = [z];   % Mean Degree

% initial conditions
s0 = 1-0.02;
b0 = 0.001;
x10 = (1-s0-b0)/2;
x20 = (1-s0-b0)/2;
rho0 = [s0 x10 x20 b0];

% CoSIS lambda parameter
lambda = 2;

% Simulation inputs
n=4;
endtime = 10;

% ---------------
% EXPERIMENTS
% ---------------

% 1: Supercritical beta

beta = 0.5/z;
DynamicsParams = [beta lambda];

% AME
[TAME xAME] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'AME');
iAME = 1-xAME(:,1);

% PA
[TPA xPA] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'PA');
iPA = 1-xPA(:,1);

% MF
[TMF xMF] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'MF');
iMF = 1-xMF(:,1);

% plot
figure
plot(TAME,iAME,'r-',TPA,iPA,'b-.',TMF,iMF,'g--','LineWidth',1)
xlabel('$t$','Interpreter','Latex');
ylabel('$i(t)\;\;\;$  ','Interpreter','Latex')
ylim([0 0.025])
l = legend('AME','PA','MF');
set(l,'Interpreter','latex','Location','northeast')
legend boxoff
set(get(gca,'ylabel'),'rotation',0)

% 1: Subcritical beta

beta = 2.0/z;
DynamicsParams = [beta lambda];

% AME
[TAME xAME] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'AME');
iAME = 1-xAME(:,1);

% PA
[TPA xPA] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'PA');
iPA = 1-xPA(:,1);

% MF
[TMF xMF] = multi_state_solver(n,DegreeDistribution,z,DynamicsParams, rho0, endtime, 'MF');
iMF = 1-xMF(:,1);

% plot
figure
plot(TAME,iAME,'r-',TPA,iPA,'b-.',TMF,iMF,'g--','LineWidth',1)
xlabel('$t$','Interpreter','Latex');
ylabel('$i(t)\;\;\;$  ','Interpreter','Latex')
ylim([0 1])
l = legend('AME','PA','MF');
set(l,'Interpreter','latex','Location','southeast')
legend boxoff
set(get(gca,'ylabel'),'rotation',0)
