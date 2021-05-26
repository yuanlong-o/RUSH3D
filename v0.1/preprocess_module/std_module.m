function [raw_stack_after_std] = std_module(preprocess_param,raw_stack,std_data_name)
%% STD Module Calculate the standard deviation of raw light field data stack with removing background
% This program is used for calculating the std value of on entire video of
% each scanning position, it also remove the background use rank 1
% normalized funciton to augment the data.

% Last update: 05/15/2021. MW

% parser
Nshift = preprocess_param.Nshift;
maxIter = preprocess_param.maxIter;

raw_size_x = size(raw_stack,1);
raw_size_y = size(raw_stack,2);

raw_stack_after_std = zeros(raw_size_x,raw_size_y,Nshift^2);

% for different scanning shift

for i = 1:Nshift^2
    tmp_stack = raw_stack(:,:,i:Nshift^2:end);  
    % get rank-1 background
    [bg_spatial, bg_temporal] = rank_1_NMF(reshape(tmp_stack,[raw_size_x*raw_size_y,size(tmp_stack,3)]), maxIter); % calculate rank 1 normalized function decomposition
    
    % calculate std
    raw_stack_after_std(:,:,i) = reshape(compute_std_image(reshape(tmp_stack, [], size(tmp_stack, 3)), bg_spatial, bg_temporal), ...
                                         size(tmp_stack, 1), size(tmp_stack, 2), []);
    disp([num2str(i), ' std has been calculated!']);
end

% imwriteTFSK(single(raw_stack_after_std),strcat(std_data_name,'.tiff'));
saveastiff(im2uint16(raw_stack_after_std / max(raw_stack_after_std(:))), strcat(std_data_name,'.tiff'))
end