function [al_pk al al_plus1 pk qk k_pk a_min a_max mult_coeff n_a_combs] = get_constants_vector(n, sys_dims, sys_dims_single_state, a_dim, k_dim, Kmin, Kmax, pkdash, combs, non_null_lin_indices, system_length)

% al_tens and al_plus1_mminus1_tens
% dimension = sys_dims = [k_dim a_dim*ones(1,n) n]
% al_tens(k,a1,..,an,l) = al if al <= node_degree, 0 otherwise
% al_plus1_tens(k,a1,..,an,l) = al+1 if al+1 <= node_degree, 0 otherwise

% al_pk_tens
% al_pk_tens(k,a1,...,an,l) = pkdask(k)*al = pk*al

al_tens = zeros(sys_dims);
al_plus1_tens = zeros(sys_dims);
al_pk_tens = zeros(sys_dims);
pk_tens = zeros(sys_dims_single_state);
a_min = zeros(k_dim,1);
a_max = zeros(k_dim,1);


for k=Kmin:Kmax
    
    k_comb = combs{k-Kmin+1};
    k_comb_indices = num2cell(k_comb+1);
    l_k_comb = size(k_comb,1);
    
    if(k == Kmin)
        a_min(k-Kmin+1) = 1;
        a_max(k-Kmin+1) = length(k_comb);
    else
        a_min(k-Kmin+1) = a_max(k-Kmin+1-1) + 1;
        a_max(k-Kmin+1) = a_max(k-Kmin+1-1) + length(k_comb);
    end
    
    for ia = 1:l_k_comb
        %al_tens(k,a1,..,an,l) = al (<=k)
        al_tens(k_comb_indices{ia,:},k-Kmin+1,:) = k_comb(ia,:);
        %al_pk_tens(k,a1,..,an,l) = al*pk (<=k)
        al_pk_tens(k_comb_indices{ia,:},k-Kmin+1,:) = k_comb(ia,:)*pkdash(k-Kmin+1);
        % al_plus1_tens(k,a1,..,an) = al+1 (if al+1 <=k)
        al_plus1_tens(k_comb_indices{ia,:},k-Kmin+1,:) = (k_comb(ia,:)+1).*(k_comb(ia,:)+1 <= k);

        
        pk_tens(k_comb_indices{ia,:},k-Kmin+1) = pkdash(k-Kmin+1);

    end

end

% % al for MF and PA, has to be sorted by degree classes
% 
% al_tens_PA_MF = permute(al_tens, [2:(n+1) 1 n+2]);
% al_tens_PA_MF = reshape(al_tens_PA_MF, [], n);
% al_PA_MF = al_tens_PA_MF(non_null_lin_indices_PA_MF, :);

% reshape the tensors into vectors

al_tens = reshape(al_tens, [], n);
al_plus1_tens = reshape(al_plus1_tens, [], n);
al_pk_tens = reshape(al_pk_tens, [], n);
pk_tens = reshape(pk_tens, [], 1);

% decompose vectors into non-null parts

al = zeros(length(non_null_lin_indices), n);
al_pk = zeros(length(non_null_lin_indices), n);
al_plus1 = zeros(length(non_null_lin_indices), n);
pk = zeros(length(non_null_lin_indices), 1);

al = al_tens(non_null_lin_indices, :);
al_pk = al_pk_tens(non_null_lin_indices, :);
al_plus1 = al_plus1_tens(non_null_lin_indices, :);
pk = pk_tens(non_null_lin_indices);


mult_coeff = zeros(system_length,1); 
% mult_coeff(ia) = k!/a_1!...a_n!. We compute this using logarithms as the
% factorials can be prohibitary large for large k
for ia=1:system_length
    k = sum(al(ia,:));
    
    mult_coeff(ia) = sum(log(1:k));
    for i=1:n
        mult_coeff(ia) = mult_coeff(ia) - sum(log(1:al(ia,i)));
    end
end
mult_coeff = exp(mult_coeff);


k_pk = zeros([n k_dim]);
for k=1:k_dim
    k_pk(:,k) = (Kmin+k-1)*pkdash(k);
end

n_a_combs = length(non_null_lin_indices);

% qk, excess degree distribution


qk = zeros(k_dim,1);
if(k_dim == 1)
    qk = 1
else
    z = dot(Kmin:Kmax,pkdash);
    for k=1:k_dim
        qk(k) = (Kmin+k-1)*pkdash(k)/z;
    end
end
    
