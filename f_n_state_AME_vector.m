function Frhs = f_n_state_AME_vector(t, x, Fvec, Rvec, al_pk, al, al_plus1, k_dim, a_dim, n, non_null_lin_indices, index_shift, index_removal)
% Version2: Instead of working with a state tensor have a vector.

x = reshape(x,[],n);

% NODE TRANSITIONS

node_transitions_out = Rvec.*x;

x_rep = repmat(x, [1 1 n]);
node_transitions_in = squeeze(sum(Fvec.*x_rep, 2));


% NEIGHBOUR TRANSITIONS

% beta_denom(l,i) is the number of links of type i-l
beta_denom = zeros(n,n);
% beta_num(l,m,i)*dt is the number of links of type i-l that change to i-m in
% an infintesmally small time step dt
beta_num = zeros(n,n,n);
beta = zeros(n,n,n);
for i=1:n
    for l=1:n
        beta_denom(l,i) = squeeze(sum(x(:,l).*al_pk(:,i)));
        % If there are no links of type i-l, stop iteration
        if(beta_denom(l,i) == 0)
            continue;
        else
            for m=1:n
                if(m == l)
                    continue;
                else
                    beta_num(l,m,i) = squeeze(sum(Fvec(:,l,m).*x(:,l).*al_pk(:,i)));
                    beta(l,m,i) = beta_num(l,m,i)/beta_denom(l,i);
                end
            end
        end
    end
end

beta_l_tot = squeeze(sum(beta, 2));
% nodes that leave the class because their neighbour changes state:
neighbour_transitions_out_rate = zeros(length(non_null_lin_indices),n);
for i=1:n
    for l=1:n
        neighbour_transitions_out_rate(:,i) = neighbour_transitions_out_rate(:,i) + al(:,l)*beta_l_tot(l,i);
    end
end
neighbour_transitions_out = neighbour_transitions_out_rate.*x;


% nodes that enter the class because their neighbour changes state:
neighbour_transitions_in = zeros(length(non_null_lin_indices),n);
for l=1:n
    for m=1:n       
        if(m==l)
            continue
        end
        
        % al -> al+1 for al<k and am -> m-1 for m>0
        x_lshift_mshift = x(index_shift{l,m},:);
        % x_lshift_mshift(al = k, am = 0) = 0
        x_lshift_mshift(index_removal{l,m},:) = 0;
        
        for i=1:n            
            neighbour_transitions_in(:,i) = neighbour_transitions_in(:,i) + beta(l,m,i)*(al_plus1(:,l).*x_lshift_mshift(:,i));
        end

    end
    
end


xdot = -node_transitions_out + node_transitions_in - neighbour_transitions_out + neighbour_transitions_in;

Frhs = reshape(xdot,[],1);

