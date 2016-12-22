function [Fvec Rvec] = get_transition_rate_vector(n, sys_dims, a_dim, k_dim, Kmin, Kmax, combs, non_null_lin_indices, DynamicsParams)

% Ftens(k,a_1,..a_n, n_1, n_2) is the rate that a k-degree nodes in state
% n_1 with a_1, ..., a_n neighbours in state 1, .., n changes to state n_2.
F_dims = [sys_dims n];
Ftens = zeros(F_dims);

for k = Kmin:Kmax
    
    k_comb = combs{k-Kmin+1};
    l_k_comb = size(k_comb,1);
    k_comb_indices = num2cell(k_comb+1);  % convert the elements of A to indices to access F
    
    
    for ia = 1:l_k_comb
        Ftens(k_comb_indices{ia,:},k-Kmin+1,:,:) = F_rates(k_comb(ia,:), DynamicsParams);
    end

    
    
end

Ftens = reshape(Ftens, [], n, n);
Fvec = zeros(length(non_null_lin_indices), n, n);
for i=1:n
    for j=1:n
        Fvec(:,i,j) = Ftens(non_null_lin_indices,i,j);
    end
end

% Rvec(k,a_1,..a_n, n) is the rate that a k-degree nodes in state
% n with a_1, ..., a_n neighbours in state 1, .., n changes to another
% state. Rtens is calculated directly from Ftens
Rvec = sum(Fvec, 3);