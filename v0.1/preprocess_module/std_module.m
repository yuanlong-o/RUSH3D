function [raw_stack_after_std] = std_module(preprocess_param,raw_stack,std_data_name)
%% STD Module Calculate the standard deviation of raw light field data stack with removing background
% This program is used for calculating the std value of on entire video of
% each scanning position, it also remove the background use rank 1
% normalized funciton to augment the data.

% Last update: 05/15/2021. MW

Nshift = preprocess_param.Nshift;


raw_size_x = size(raw_stack,1);
raw_size_y = size(raw_stack,2);

mean_data = zeros(raw_size_x,raw_size_y,Nshift^2);

var_data = zeros(raw_size_x,raw_size_y,Nshift^2);
raw_stack_after_std = zeros(raw_size_x,raw_size_y,Nshift^2);

for i = 1:Nshift^2
    tmp_stack = raw_stack(:,:,i:Nshift^2:end);  
    [bg_spatial, bg_temporal] = rank_1_NMF(reshape(tmp_stack,[raw_size_x*raw_size_y,size(tmp_stack,3)]), maxIter); % calculate rank 1 normalized function decomposition
    tmp_stack_tmp = tmp_stack - reshape(bg_spatial*bg_temporal,[raw_size_x,raw_size_y,size(tmp_stack,3)]); 
    tmp_stack = abs(tmp_stack -  tmp_stack_tmp); % detrending
    clear bg_spatial  bg_temporal;
    tmp_stack(tmp_stack<0) = 0;
    
    mean_data(:,:,i) = sum(tmp_stack,3)./size(tmp_stack,3);
    var_data(:,:,i) =  sum((tmp_stack- repmat(mean_data(:,:,i),[1,1,size(tmp_stack,3)])).^2,3);
    raw_stack_after_std(:,:,i) = sqrt(double(var_data(:,:,i))./(size(tmp_stack,3)-1)); % calculate variance and std
    disp([num2str(i), ' std has been calculated!']);
end

imwriteTFSK(single(raw_stack_after_std),strcat(std_data_name,'.tiff'));
end