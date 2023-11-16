function spatial_comp_render(A, A_manual, dims)
curr_comp_crop_array = [];
for i = 1 : size(A, 2)
    curr_comp = full(A(:, i));
    curr_comp = reshape(curr_comp, dims(1), dims(2));
    curr_mask = reshape(A_manual(:, i), dims(1), dims(2));

    
    top_left = min(find(curr_mask ));
    bottom_right = max(find(curr_mask ));
    
%     center = (top_left + bottom_right) / 2;
    [top_left_h, top_left_w] = ind2sub([dims(1), dims(2)], top_left);
    [bottom_right_h, bottom_right_w] = ind2sub([dims(1), dims(2)], bottom_right);
%     top_left_h = top_left_h - round(gSiz / 2);
%     top_left_w = top_left_w - round(gSiz / 2);
%     bottom_right_h = top_left_h + gSiz;
%     bottom_right_w = top_left_w + gSiz;
    
    curr_comp_crop = curr_comp(top_left_h : bottom_right_h, top_left_w : bottom_right_w);
    curr_comp_crop = curr_comp_crop / max(curr_comp_crop(:));
    
    curr_comp_crop_array{i} = curr_comp_crop; 
end
montage(curr_comp_crop_array)
end