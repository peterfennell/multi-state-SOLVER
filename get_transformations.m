function [combs non_null_lin_indices full_to_reduced_map system_length] = get_transformations(k_dim, a_dim, n, Kmin, Kmax, sys_dims_single_state)

non_null_tens = zeros(sys_dims_single_state);

% Now, for each possible value of the degree k, calculate all possible
% combinations a1,...,an such that a1+...+an=k. Mark these elements with a
% one. All the other elements are null-elements

% First, calculate a cell such that for each value of k, we get all
% possible combinations of n elements that sum to k
combs = cell(k_dim);

for k=Kmin:Kmax

    % This next part caculates the all possible combinations a1,..,an that sum
    % to k
 
    % explanation: divide the k units in the sum into n intervals,
    % seperated by n-1 dividers. Then a1 is the number between the start and the 1st divider, ai is the number between the i-1st divider and the ith divider, etc.  Get all such combinations:
    c = nchoosek(1:(n+k-1),n-1);         % All combinations of n-1 elements from the vector [1 2 3 ... n+k-1]
    c_size = size(c,1);                  % The number of them combinations
    A = zeros(c_size,n);                 % The matrix of all possible combinations a_1,...,a_n such that \sum_i a_i = k
    % add divider at the start and at the end:
    c = cat(2, zeros(c_size, 1), c, (n+k)*ones(c_size, 1));
    % subtracting dividers gives the number of elements between each divider
    k_comb = diff(c, 1, 2) - 1;
    
    combs{k-Kmin+1} = k_comb;
    
end
    
% Next, for each of those elements mark a one in the
% non_null_elements_matrix

for k=Kmin:Kmax

    k_comb = combs{k-Kmin+1};
    k_comb_indices = num2cell(k_comb+1);  % convert the elements of k_comb to indices to access F

    
    %non_null_tens(k,a1,..,an) = 1 if the combination of coefficients is
    %possible i.e. a1+...+an=k
    for ia=1:length(k_comb)
        non_null_tens(k_comb_indices{ia,:},k-Kmin+1) = 1;
    end

end

% Now, compute the linear indices of the non-zero terms
non_null_lin_indices = find(non_null_tens);
system_length = length(non_null_lin_indices);
% non_null_lin_indices(i) is the linear index of the tensor element. on the
% other hand, need to calculate i if we know the linear index
full_to_reduced_map = zeros(k_dim*(a_dim^n),1);
full_to_reduced_map(non_null_lin_indices) = 1:system_length;

% % Same for PA and MF quantities, need them sorted by degree, though
% 
% tensor_index_dims_PA_MF = [a_dim*ones(1,n) k_dim];
% non_null_tens_PA_MF = zeros(tensor_index_dims_PA_MF);
% 
% for k=Kmin:Kmax
% 
%     k_comb = combs{k-Kmin+1};
%     k_comb_indices = num2cell(k_comb+1);  % convert the elements of k_comb to indices to access F
% 
%     
%     %non_null_tens(k,a1,..,an) = 1 if the combination of coefficients is
%     %possible i.e. a1+...+an=k
%     for ia=1:length(k_comb)
%         non_null_tens_PA_MF(k_comb_indices{ia,:},k-Kmin+1) = 1;
%     end
% 
% end
% 
% % Now, compute the linear indices of the non-zero terms
% non_null_lin_indices_PA_MF = find(non_null_tens_PA_MF);



