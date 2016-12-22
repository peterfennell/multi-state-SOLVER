function [t_points x_tots] = multi_state_solver(n, DegreeDistribution,DistParams,DynamicsParams,rho0,endtime, scheme)

%-----------------------------------------------------------
%% INPUTS
% n: number of states
% Degree distribution:
% - choices:'zRegular'  (degree-regular/Bethe lattice)
%           'PRG'       (poisson random graph (ER network))
%           'truncSFN'  (truncated scale-free network)
%           'custom'
% DistParams: parameters of the Degree distrubition. Here z = mean degree, Kmin = min degree, Kmax = max degree.
% -         if 'zRegular': DistParams = [z];
%           if 'PRG'     : DistParams = [Kmin Kmax z];
%           if 'truncSFN': DistParams = [Kmin Kmax gamma]; (gamma = power law exponent)
%           if 'custom'  : DistParams(i) = Pr(degree=i) for 0 <= i <= Kmax (e.g., DistParams = [0 0 0 0 1] for Degree-regular network with z=4)
% DynamicsParams: optional vector of parameters to be passed to the F_rates.m function
% rho0: vector of initial conditions where rho0(i) = inital fraction of nodes in the network in state i.
% -         rho0 must satisfy length(rho0) = n and sum(rho0)=1;
% endtime: system is soved for times 0 -> endtime
% scheme: Which approximation scheme to employ
% - choices:'MF' (Mean-field)
%           'PA' (Pair approximatio)
%           'AM' (Approximate Master Equation)
%
%-----------------------------------------------------------
%% OUTPUT
% t_points: Time points t_j for which solution is given for 1 <= j <= 100, 0 <= t_j <= endtime
% x_tots:   x_tots(t_j,i) is the fraction of nodes in the network in state i at time t_j
%
%-----------------------------------------------------------

%% PROGRAM

% Get the degree distribution Pk aswell as Kmin, Kmaz, average degree etc. DegreeDistribution = 'zRegular'; % 'PRG' or 'truncSFN' or 'zRegular' or 'custom'

[Kmin, Kmax, z, gamma, pkdash] = get_degree_distribution(DegreeDistribution,DistParams);


% Various Dimensions
k_dim = Kmax-Kmin+1;        % number of possible values for k
a_dim = Kmax+1;             % number of possible values for the a's (0<=a<=k)

% dimensions of the system i.e. all possible combinations of (k,a_1,...a_n)
sys_dims = [repmat(a_dim,1,n) k_dim n];   
sys_dims_single_state = [repmat(a_dim,1,n) k_dim];

% The system we are studying can be represented as a tensor, where the
% variables x(k,a1,..,an,n) are the fraction of k-degree nodes in the
% network in state n with a1,..,an neighbours in states 1,..,n
% respectively. However, this is subjust to the constraints that
% a1+...+an=k and any combination of ai's that doesn't satisfy this will
% just have a null entry in the tensor. This is extremely wasteful both in
% terms of memory and speed.

% To combat this we reduce the system. Firstly, we transform the tensor
% into a 2-d system x(j(k,a1,..,an),n), where j(k,a1,..,an) is the linear
% index of the tensor index [k a1 ... an]. This j(k,a1,..,an) index still
% contains the null elements of system, and so we further secompose 
% j(k,a1,..,an) into j(i), where j(i) = j(k_i,a1_i,..,an_i) corresponding
% to the non-null elements of the system. The first part of the programme
% is to calculate these transformations

[combs non_null_lin_indices full_to_reduced_map system_length] = get_transformations(k_dim, a_dim, n, Kmin, Kmax, sys_dims_single_state);


% Calculate the constants that are needed in the AME framework. These
% include al, al_plus1 and al_pk which signify al(k,a1,..,an,l) = al,
% al_plus1(k,a1,..,an,l) = al+1 and al_pk(k,a1,..,an,l) = pkdash(k)*al
% respectively

[al_pk al al_plus1 pk_AME qk k_pk a_min a_max mult_coeff n_a_combs] = get_constants_vector(n, sys_dims, sys_dims_single_state, a_dim, k_dim, Kmin, Kmax, pkdash, combs, non_null_lin_indices, system_length);


% Calculate the transofmation vectors that will transform from
% (k,a1,...,an) to (k,a1,..,al+1,..,am-1,..,an)

[index_removal index_shift] = get_shifting_vectors(Kmin, Kmax, k_dim, a_dim, n, combs, non_null_lin_indices, full_to_reduced_map, sys_dims_single_state);

% state tensor: x(k,a_1,..,a_n,n) is the fraction of k-degree nodes in state
% n with a_1, ..., a_n neighbours in state 1, .., n respectively 
% Get the initial fraction of nodes in each compartment
if(sum(rho0) ~= 1)
    disp(strcat('Sum of rho0 is NOT one; = ',num2str(sum(rho0))));
    s=-1;
    return;
end

[x0_AME x0_PA x0_MF] = initial_conditions_vector(rho0, n, sys_dims, a_dim, a_min, k_dim, Kmin, Kmax, combs, non_null_lin_indices, mult_coeff);

% Transition rate tensor: Ftens(k,n_1,n_2,a_1,..a_n) is the rate that a k-degree nodes in state
% n_1 with a_1, ..., a_n neighbours in state 1, .., n changes to state n_2.
% Rtens(k,n,a_1,..a_n) is the rate that a k-degree nodes in state
% n with a_1, ..., a_n neighbours in state 1, .., n changes to another
% state. Rtens is calculates directly from Ftens
[Fvec Rvec] = get_transition_rate_vector(n, sys_dims, a_dim, k_dim, Kmin, Kmax, combs, non_null_lin_indices, DynamicsParams);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ODE SOLVER
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set ODE solver parameters
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1.5); % defaults are 'RelTol',1e-3,'AbsTol',1e-6
tspan= [0 endtime]; %0:endtime/100:endtime;
% tspan = logspace(-2,log10(endtime),300); %0:endtime/100:endtime; %


% SOLVERS



% MF

if(scheme == 'MF')
    fprintf('System initialized, MF simulation commenced.\n');
    [T_MF,Y_MF]=ode2r(@(t,x) f_n_state_MF_vector(t, x, Fvec, Rvec, al, k_dim, Kmin, n, a_min, a_max, k_pk, qk, mult_coeff, n_a_combs), tspan, x0_MF, options);
    Y_MF = reshape(Y_MF,[length(T_MF) k_dim n]);
    X_MF = zeros(length(T_MF),n);
    for i=1:n
        for t=1:length(T_MF)
            if(k_dim > 1)
                X_MF(t,i) = dot(Y_MF(t,:,i),pkdash);
            else
                X_MF(t,i) = Y_MF(t,:,i);
            end
        end
    end
    fprintf('Numerical solution completed\n');
    t_points = T_MF;
    x_tots   = X_MF;
elseif(scheme == 'PA')
    fprintf('System initialized, PA simulation commenced.\n');
    [T_PA,Y_PA]=ode2r(@(t,x) f_n_state_PA_vector(t, x, Fvec, Rvec, al_pk, al, k_dim, Kmin, n, a_min, a_max, k_pk, mult_coeff, n_a_combs), tspan, x0_PA, options);
    X_PA_temp = reshape(Y_PA(:,1:(n*k_dim)),[length(T_PA) n k_dim]);
    X_PA = zeros(length(T_PA),n);
    for i=1:n
        for t=1:length(T_PA)
            if(k_dim > 1)
                X_PA(t,i) = dot(squeeze(X_PA_temp(t,i,:)),pkdash);
            else
                X_PA(t,i) = squeeze(X_PA_temp(t,i,:));
            end
        end
    end
    Q_PA = reshape(Y_PA(:,(n*k_dim+1):end),[length(T_PA) n n k_dim]);
    fprintf('Numerical solution completed\n');
    t_points = T_PA;
    x_tots   = X_PA;
else
    fprintf('System initialized, AME solver commenced.\n');
    [T_AME,Y_AME]=ode2r(@(t,x) f_n_state_AME_vector(t, x, Fvec, Rvec, al_pk, al, al_plus1, k_dim, a_dim, n, non_null_lin_indices, index_shift, index_removal), tspan, x0_AME, options);
    Y_AME  = reshape(Y_AME, [length(T_AME) length(non_null_lin_indices) n]);
    x_tots = zeros(length(T_AME), n);
    for t=1:length(T_AME)
        for i=1:n
            X_AME(t,i) = dot(Y_AME(t,:,i),pk_AME);
        end
    end
    fprintf('Numerical solution completed\n');
    t_points = T_AME;
    x_tots   = X_AME;
end


