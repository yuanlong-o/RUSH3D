function patch_info_array = determine_patch_size(size_h, size_w, estimate_patch_size)
    M = ceil(size_h / estimate_patch_size);
    N = ceil(size_w / estimate_patch_size);
    
    refined_patch_size_h = ceil(size_h / M);
    refined_patch_size_w = ceil(size_w / N);
    
    % genrate patch size
    patch_info_array = [];
    for i = 1 : M 
        for j = 1 : N
            top_left = [(i - 1) * refined_patch_size_h + 1, (j - 1) * refined_patch_size_w + 1];

            % patch location
            % due to the ceil option in the above, this is safe
            bottom_right = [min(i * refined_patch_size_h, size_h), ...
                min(j * refined_patch_size_w, size_w)];
            
            % udpate patch size
            patch_size_h_new = bottom_right(1) - top_left(1) + 1;
            patch_size_w_new = bottom_right(2) - top_left(2) + 1;

            patch_loc = [top_left; bottom_right];      

            
            patch_info.location = patch_loc;
            patch_info.size = [patch_size_h_new, patch_size_w_new];
            
            patch_info_array_length = length(patch_info_array);
            patch_info_array{patch_info_array_length + 1} = patch_info;
        end
    end
end