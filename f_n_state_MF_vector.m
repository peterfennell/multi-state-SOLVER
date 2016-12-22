function f_rhs = f_n_state_MF_vector(t, x, Fvec, Rvec, al, k_dim, Kmin, n, a_min, a_max, k_pk, qk, mult_coeff, n_a_combs)
% Version2: Instead of working with a state tensor have a vector.

% Variables: 
% x_MF(i,k): expected fraction of k-degree nodes in the network that are in
% state i 

x_MF = reshape(x,[k_dim n]);


% STEP 1
% Compute Mult(ia), the multinomial probability that a k(=|a(ia)|)-degree node in state i
% has a(ia) neighbours in the various different states

Mult = zeros(n_a_combs);
omega = zeros(n,1); % omega(j) is the probability that the neighbour of a node is in state j
for j=1:n
    omega(j) = dot(qk,x_MF(:,j));
end
omega = transpose(repmat(omega, 1, n_a_combs));

Mult = prod(squeeze(omega.^al),2).*mult_coeff;
Mult = repmat(Mult,1,n);

F_times_x_out = zeros(n_a_combs,n);
F_times_x_in = zeros(n_a_combs,n);
% F_times_x_MF(ai,i,j) = Fvec(ai,i,j)*x_MF(i,k);


for k=1:k_dim
    
    for i=1:n
        
        F_times_x_out(a_min(k):a_max(k),i) = squeeze(sum(Fvec(a_min(k):a_max(k),i,:),3)).*x_MF(k,i);
        
        for j=1:n
            F_times_x_in(a_min(k):a_max(k),i) = F_times_x_in(a_min(k):a_max(k),i) ...
                + Fvec(a_min(k):a_max(k),j,i).*x_MF(k,j);
        end
    end
end

f_rhs_presum = zeros(n_a_combs,n);
f_rhs_presum = -(F_times_x_out - F_times_x_in).*Mult;

f_rhs = zeros(k,n);
for k=1:k_dim
    f_rhs(k,:) = sum(f_rhs_presum(a_min(k):a_max(k),:),1);
end

f_rhs = reshape(f_rhs,[],1);

