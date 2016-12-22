function [x0_AME x0_PA x0_MF] = initial_conditions_vector(rho0, n, sys_dims, a_dim, a_min, k_dim, Kmin, Kmax, combs, non_null_lin_indices, mult_coeff)
% at some stage introduce fix to calculate reduce problem in the case that
% one of the initial conditions is 0

x0_tens = zeros(sys_dims);



for k=Kmin:Kmax

    k_comb = combs{k-Kmin+1};
    l_k_comb = size(k_comb,1);
    k_comb_indices = num2cell(k_comb+1);

    a_min_k = a_min(k-Kmin+1);
    
    for ia = 1:l_k_comb
        % Multinomial probability of having a_i neighbours in the
        % state i
        x0_tens(k_comb_indices{ia,:},k-Kmin+1,:) = rho0*mult_coeff(a_min_k+ia-1)*(prod(rho0.^k_comb(ia,:)));
    end

    
end

x0_tens = reshape(x0_tens,[],n);
x0_AME(:,:) = x0_tens(non_null_lin_indices,:);
x0_AME = reshape(x0_AME, [], 1);


x0_PA = zeros([n k_dim]);
q0_PA = zeros([n n k_dim]);

for i=1:n
    x0_PA(i,:) = rho0(i);
    for j=1:n
        q0_PA(j,i,:) = rho0(j);
    end
end

x0_PA = reshape(x0_PA,[],1);
q0_PA = reshape(q0_PA,[],1);

x0_PA = cat(1,x0_PA, q0_PA);


x0_MF = zeros([k_dim n]);

for i=1:n
    x0_MF(:,i) = rho0(i);
end

x0_MF = reshape(x0_MF, [], 1);

