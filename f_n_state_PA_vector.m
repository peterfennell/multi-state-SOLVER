function f_rhs = f_n_state_PA_vector(t, x, Fvec, Rvec, al_pk, al, k_dim, Kmin, n, a_min, a_max, k_pk, mult_coeff, n_a_combs)
% Version2: Instead of working with a state tensor have a vector.

% Variables: 
% x_PA(i,k): expected fraction of k-degree nodes in the network that are in
% state i 
% q_PA(j,i,k): expected fraction of state j nodes at the end of links
% emenating from k-degree nodes in state i

x_PA = reshape(x(1:(n*k_dim)),[n k_dim]);
q_PA = reshape(x((n*k_dim+1):end,1),[n n k_dim]);

% STEP 1
% Compute Mult_times_x(a,i) = Mult(a,i)*x_PA(i,k(=sum(a))), where Mult(ia, i)
% is the multinomial probability that a k(=|a(ia)|)-degree node in state i
% has a(ia) neighbours in the various different states

Mult_times_x = zeros(n_a_combs, n);
q_rep = repmat(q_PA,[1 1 1 n_a_combs]);
q_rep = permute(q_rep, [4 1 2 3]);

for k=1:k_dim
    
    for i=1:n      
        Mult_times_x(a_min(k):a_max(k),i) =  x_PA(i,k)* ...
            prod(squeeze(q_rep(a_min(k):a_max(k),:,i,k)).^al(a_min(k):a_max(k),:),2) ...
            .*mult_coeff(a_min(k):a_max(k));
    end
       
end

clear('q_rep');    % Clear up some memory



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NODE EQUATION
% Evolution equation for x_PA(i,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_rhs = zeros(n,k_dim);

x_out_presum = Rvec.*Mult_times_x;
x_in_presum = zeros(n_a_combs, n);
for i=1:n
    for j=1:n
        if(j==i)
            continue
        else
            x_in_presum(:,i) = x_in_presum(:,i) + Fvec(:,j,i).*Mult_times_x(:,j);
        end
    end
end

x_rhs_presum = - x_out_presum + x_in_presum;
for k=1:k_dim
    x_rhs(:,k) = sum(x_rhs_presum(a_min(k):a_max(k),:),1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINK EQUATION
% Evolution equation for x_PA(i,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_rhs = zeros(n,n,k_dim);


% beta_denom(i,l) is the number of links of type i-l
beta_denom = zeros(n,n);
% beta_num(i,l,m)*dt is the expected number of links of type i-l that change to i-m in
% an infintesmally small time step dt
beta_num = zeros(n,n,n);
beta = zeros(n,n,n);
for i=1:n
    
       
    for l=1:n
    
        beta_denom(i,l) = sum(x_PA(l,:).*k_pk(l,:).*transpose(squeeze(q_PA(i,l,:)))); %k_pk(l,k)
        
        % If there are no links of type i-l, stop iteration
        if(beta_denom(i,l) == 0)
            continue;
        else
            for m=1:n
                if(m == l)
                    continue;
                else
                    beta_num(i,l,m) = squeeze(sum(Fvec(:,l,m).*Mult_times_x(:,l).*al_pk(:,i)));
                    beta(i,l,m) = beta_num(i,l,m)/beta_denom(i,l);
                end
            end
        end
    end
end


beta_out_tot = transpose(squeeze(sum(beta, 3)));
beta_permuted = permute(beta,[2 1 3]);

% beta_out_tot(j,i) is rate neighbour of a node in state i changes from state
% j



q_rhs = zeros(n,n,k_dim);
q_neighb_out = zeros(n,n,k_dim);
q_neighb_in = zeros(n,n,k_dim);

for k=1:k_dim
    q_neighb_out(:,:,k) = q_PA(:,:,k).*beta_out_tot;
    for l=1:n
        for m=1:n
            q_neighb_in(l,:,k) = q_neighb_in(l,:,k) + q_PA(m,:,k).*beta_permuted(m,:,l);
        end
    end
end

q_node_change = zeros(n,n,k_dim);

for k=1:k_dim
    for i=1:n
        if(x_PA(i,k) == 0)
            continue;
        end
        for j=1:n
            q_node_change(j,i,k) = - (1/x_PA(i,k))* ...
                sum(x_rhs_presum(a_min(k):a_max(k),i).*(q_PA(j,i,k) - (al(a_min(k):a_max(k),j)/(Kmin+k-1))));
        end
    end
end

q_rhs = q_node_change - q_neighb_out + q_neighb_in;



x_rhs = reshape(x_rhs,[],1);
q_rhs = reshape(q_rhs,[],1);


f_rhs = [x_rhs; q_rhs];

