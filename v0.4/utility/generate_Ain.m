function A_in = generate_Ain(seg, center, size_h, size_w)
    A_in = zeros(size_h, size_w, length(seg));
    for i = 1 : length(seg)
        curr_center = center(i, 1 : 2);
        [curr_patch_size_h, curr_patch_size_w] = size(seg{i});
        curr_patch_top_left = [curr_center(1) - floor(curr_patch_size_h / 2), ...
                               curr_center(2) - floor(curr_patch_size_w / 2)];
        curr_patch_top_left(curr_patch_top_left < 1) = 1;

        curr_patch_bottom_right = [curr_patch_top_left(1) + curr_patch_size_h - 1, ...
                               curr_patch_top_left(2) + curr_patch_size_w - 1];
        
        if curr_patch_bottom_right(1) > size_h;  curr_patch_bottom_right(1) = size_h; end
        if curr_patch_bottom_right(2) > size_w;  curr_patch_bottom_right(2) = size_w; end
            
        update_patch_size_h = curr_patch_bottom_right(1) - curr_patch_top_left(1) + 1;
        update_patch_size_w = curr_patch_bottom_right(2) - curr_patch_top_left(2) + 1;
        
        A_in(curr_patch_top_left(1) : curr_patch_bottom_right(1),...
             curr_patch_top_left(2) : curr_patch_bottom_right(2),...
             i) = seg{i}(1 : update_patch_size_h, 1 : update_patch_size_w);
    end
    % average the patch
    A_in = mean(A_in, 3);
end
