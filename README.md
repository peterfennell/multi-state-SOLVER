# multi-state-SOLVER
MATLAB/Octave code for numerically solving multi-state degree-based frameworks, as presented in [1]. 

## Program Files

The main function is multi_state_solver(). This has the form

**function [t_points x_tots] = multi_state_solver(n,DegreeDistribution,DistParams,DynamicsParams,rho0,endtime,scheme)**,

and the various function inputs and outputs are described in lines 3-30 of the file multi_state_solver.m. 

The rate matrix function F for the dynamics (see [1]) should be specified in the file F_rates.m. Example rate functions for different dynamics including the SI and SIS model, the Bass diffusion model, the generalize model of interacting diseases of Sanz et al. [2] and the Ising Glauber model are included. 

## Examples

Two comprehensive examples are included here; these are ExCooperativeSISmodel.m which reproduces Fig. 4 in [1] and ExFAmodel.m and which reproduces Fig. 8 in [1]. These examples should make the usage of the multi_state_solver very clear. When using these examples please make sure to uncomment the relative rate functions in the F_rates.m files (this is explained in the example .m files). 

## Notes

Note the following three important points:
 - The pair approximation ("PA") scheme requires non-zero initial conditions. The reason for this is a division by zero of the state variables (as seen in Eq.(18) of [1]).
 - The default ode solver is used here, ode45. For stiff systems, ode2r should be employed to gaurantee accuracy (This can be done by replacing "ode2r" with "ode45" in lines 116 (MF), 133 (PA) and/or 151 (AME) of the file multi_state_solver.m)
 - If alternative time stamps are required (from the output t_points), specific time stamps can be specified by the tspan argument on line 104 of multi_state_solver.m. For example, the commented-out line 105 gives time stamps that are evenly spaced, while the commented-out line 106 gives time stamps that are evenly spaced on a logaritmic x-axis. 

[1] Fennell, P.G., Gleeson, J.P., "Multistate dynamical processes on networks: Analysis through degree-based approximation frameworks.", 	arXiv:1709.09969

[2] Sanz et al. "Dynamics of Interacting Diseases",Phys Rev X, 2014
