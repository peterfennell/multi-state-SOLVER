function [index_removal index_shift] = get_shifting_vectors(Kmin, Kmax, k_dim, a_dim, n, combs, non_null_lin_indices, full_to_reduced_map, sys_dims_single_state)

% Index shifters
index_shift = cell(n,n);

% each index (k,a1,..,an) of a tensor has a corresponding decomposed linear index
% j(i). The question is what is the decomposed linear index of
% (k,a1,al+1,am-1,an) for given values of l and m 

tensor_lin_indices = (1:k_dim*(a_dim^n))';
tensor_lin_indices = reshape(tensor_lin_indices, sys_dims_single_state);
tensor_lin_indices_orig = tensor_lin_indices(non_null_lin_indices);


for l=1:n
    for m=1:n
        if(l==m)
            continue
        end
        % shift the tensor +1 in the l-direction and -1 in the m-direction
        tensor_lin_indices_lpos_shifts = circshift(tensor_lin_indices, [zeros(1,l-1) -1 zeros(1,n-l)]);
        tensor_lin_indices_lpos_mneg_shifts = circshift(tensor_lin_indices_lpos_shifts, [zeros(1,m-1) +1 zeros(1,n-m)]);
    
        % reshape the tensor into a column vector
        tensor_lin_indices_lpos_mneg_shifts = reshape(tensor_lin_indices_lpos_mneg_shifts, [], 1);
        % take only the system parts
        tensor_lin_indices_lpos_mneg_shifts_vec = tensor_lin_indices_lpos_mneg_shifts(non_null_lin_indices);
    
        % now compute what elemnts of the reduced vector have changed and
        % where they have changed to. This is incoded in index_shift_vec;
        % if index_shift_vec(i) = i then there has not been a shift,
        % otherwise it has moved
        index_shift_vec = (1:length(non_null_lin_indices))';
        for i=1:length(non_null_lin_indices)
            
            % check if the index of the shifted, decomposed vector has
            % changed
            if(tensor_lin_indices_lpos_mneg_shifts_vec(i) ~= tensor_lin_indices_orig(i))
                % if it has, then find the position that it has moved to...
                orig_index = full_to_reduced_map(tensor_lin_indices_lpos_mneg_shifts_vec(i));
                if(orig_index == 0)
                    % the shift has mapped it to a null combination of the
                    % system
                    continue
                else
                    index_shift_vec(i) = orig_index;
                end
            end
        end

        index_shift{l,m} = index_shift_vec;
        
    end
end




% Index removal
index_removal = cell(n,n);

% certain shifts are not possible. If l=k then there can't be a shift to
% l+1. if m=0 then there cant be a shift to m=-1.



for k=Kmin:Kmax
 
    k_comb = combs{k-Kmin+1};
    l_k_comb = size(k_comb,1);
    k_comb_indices = num2cell(k_comb+1);
    
    for l=1:n
        for m=1:n
            
            if(m==l)
                continue;
            end
            
            index_removal_tens = zeros(sys_dims_single_state);

            nulls = find((k_comb(:,l) == k)|(k_comb(:,m) == 0));
            l_nulls = size(nulls,1);

            for ia = 1:l_nulls
                index_removal_tens(k_comb_indices{nulls(ia),:},k-Kmin+1) = 1;
            end

            % linear_indices of elements to be removed
            index_removal_vec = find(index_removal_tens);

            % decomposed indices of elements to be removed;
            decomposed_index_removal_vec = full_to_reduced_map(index_removal_vec);

            index_removal{l,m} = [index_removal{l,m}; decomposed_index_removal_vec];


        end
    end



end
