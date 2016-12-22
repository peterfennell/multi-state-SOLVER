function [Kmin, Kmax, z, gamma, pkdash] = get_degree_distribution(DegreeDistribution,DistParams) 

z=0;
gamma=0;

if strcmpi(DegreeDistribution,'PRG')
    Kmin=DistParams(1);
    Kmax=DistParams(2);
    z=DistParams(3);
    kv=(Kmin:Kmax)'; % column vector, K-Kmin+1 by 1.
    pkdash=z.^kv./(factorial(kv))*exp(-z);
    pkdash=pkdash./sum(pkdash);
elseif strcmpi(DegreeDistribution,'zRegular')|strcmpi(DegreeDistribution,'Bethe')
    Kmin=DistParams(1);
    Kmax=DistParams(1);
    z=DistParams(1);
    kv=(Kmin:Kmax)'; % column vector, K-Kmin+1 by 1.
    pkdash=(kv==z); % z-regular (Bethe lattice)
elseif strcmpi(DegreeDistribution,'truncSFN')
    Kmin=DistParams(1);
    Kmax=DistParams(2);
    gamma=DistParams(3);
    kv=(Kmin:Kmax)'; % column vector, K-Kmin+1 by 1.
    Pk_unnorm=kv.^gamma;
    pkdash=Pk_unnorm/sum(Pk_unnorm);
elseif strcmpi(DegreeDistribution,'custom')
    Kmin   = find(DistParams,1)-1;
    pkdash = DistParams((find(DistParams,1)):end);
    Kmax   = Kmin+length(pkdash)-1;
else
    disp('Degree distribution not recognized.')
end 

if(sum(pkdash) ~= 1)
    disp(strcat('ERROR: Sum of Pk should be 1; = ',num2str(sum(pkdash))));
    disp(strcat('Mean degree= ',num2str(pkdash'*kv)));
end