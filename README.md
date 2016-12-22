# multi-state-SOLVER
MATLAB/Octave code for numerically solving multi-state degree-based frameworks, as presented in [1]. 

The main function is multi_state_solver(). This has the form

**function [t_points x_tots] = multi_state_solver(n,DegreeDistribution,DistParams,DynamicsParams,rho0,endtime,scheme)**,

and the various function inputs and outputs are described in lines 3-30 of the file multi_state_solver.m. 

The rate matrix function F for the dynamics (see [1]) should be specified in the file F_rates.m. Example rate functions for different dynamics including the SI and SIS model, the Bass diffusion model, the generalize model of interacting diseases of Sanz et al. [2] and the Ising Glauber model are included. 

For example, to reproduce Fig. 4 (left) in [1], uncomment lines 43-48 of F_rates.m and run the command

**multi_state_solver(4,'zRegular',4,[0.0625 2],[0.9779 0.01 0.01 0.001],10, scheme)**

in MATLAB or Octave, where scheme = 'MF' for the MF curve, 'PA' for the PA curve of 'AM' for the AME curve

Note the following two important points:
 - By default, an ode solver for stiff system is employed to gaurantee accuracy. This may be uneccesary for many dynamics, in which cases ode45 should be used (by replacing "ode2r" with "ode45" in lines 116 (MF), 133 (PA) and/or 155 (AME) of the file multi_state_solver.m)
 - If alternatice time stamps are required (from the output t_points), specific time stamps can be specified by the tspan argument on line 104 of multi_state_solver.m. For example, the commented-out line 105 gives time stamps that are evenly spaces on a logaritmic x-axis. 

[1] Fennell, P.G., Gleeson, J.P., "Multistate dynamical processes on networks: analysis through degree-based approximation frameworks.", in preparation

[2] Sanz et al. "Dynamics of Interacting Diseases",Phys Rev X, 2014
