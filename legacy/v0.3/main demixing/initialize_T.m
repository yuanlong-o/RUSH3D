function [T_init, T_shell_init] = initialize_T(processed_video, S_mask_init,  S_shell_mask_init)
%% this function is used to initialize T
%  last update: 5/23/2021. YZ

[size_h, size_w, valid_frame, valid_wigner] = size(processed_video);
        
valid_comp = size(S_mask_init, 1);
T_init = zeros(valid_comp, size(processed_video, 3));
T_shell_init = zeros(valid_comp, size(processed_video, 3));


for i = 1 : valid_comp
    if mod(i, 10) == 0
       fprintf('%d comp in %d \n', i, valid_comp) 
    end
    for j = 1 : valid_wigner
        curr_S = S_mask_init{i, j}; % broadcast here
        curr_shell = S_shell_mask_init{i, j};
        curr_view = processed_video(:, :, :, j);
        buf1(:, j) =squeeze(mean( bsxfun(@times, curr_view, full(curr_S)), [1, 2]));
        buf2(:, j) =squeeze(mean( bsxfun(@times, curr_view, full(curr_shell)), [1, 2]));        
        
    end
    T_init(i, :) = mean(buf1, 2); % average the signals from different wigner
    T_shell_init(i, :) = mean(buf2, 2);
end
end