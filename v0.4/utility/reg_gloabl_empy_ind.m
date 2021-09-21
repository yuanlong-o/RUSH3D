function global_reg = reg_gloabl_empy_ind(empty_ind_global, empty_ind)
% empty_ind_global : 0-1 array. can be longer
% empty_ind: 0-1 array
% global_reg : 0-1 array

pre_valid = find(empty_ind_global);
curr_valid = pre_valid(find(empty_ind == 0));
global_reg = zeros(1, length(empty_ind_global));
global_reg(curr_valid) = 1;
global_reg = 1 - global_reg;


end