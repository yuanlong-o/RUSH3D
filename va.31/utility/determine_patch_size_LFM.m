function [patch_info_array] = determine_patch_size_LFM(size_hv, size_wv, size_h_wdf, size_w_wdf, seed_param)

Nnum = seed_param.Nnum;
estimate_patch_size = seed_param.estimate_patch_size;
vol_ds = seed_param.vol_ds;
wdf_ds = seed_param.wdf_ds;
NN_ds = seed_param.NN_ds;
overlap = seed_param.overlap;
overlap_ds = overlap/wdf_ds;
ds_ori_ratio = NN_ds / vol_ds;
% ds_ratio = NN_ds * wdf_ds / vol_ds;

M = ceil(size_hv / estimate_patch_size);
N = ceil(size_wv / estimate_patch_size);

refined_patch_size_h = ceil(size_hv / M);
refined_patch_size_h = refined_patch_size_h - mod(refined_patch_size_h, Nnum) + Nnum;

refined_patch_size_w = ceil(size_wv / N);
refined_patch_size_w = refined_patch_size_w - mod(refined_patch_size_w, Nnum)+ Nnum;

% genrate patch size
patch_info_array = [];
for i = 1 : M
    for j = 1 : N
        top_left = [(i - 1) * refined_patch_size_h + 1, (j - 1) * refined_patch_size_w + 1];
        top_left_wdf = max(ceil(top_left/ds_ori_ratio),1) + overlap;
        top_left_wdf_ds = max(ceil(top_left/(ds_ori_ratio*wdf_ds)),1) + overlap_ds;
        % patch location
        % due to the ceil option in the above, this is safe
        bottom_right = [min(i * refined_patch_size_h, size_hv), ...
            min(j * refined_patch_size_w, size_wv)];
        bottom_right_wdf = min(ceil(bottom_right/ds_ori_ratio) + overlap,[size_h_wdf, size_w_wdf]);
        bottom_right_wdf_ds = min(ceil(bottom_right/(ds_ori_ratio*wdf_ds)) + overlap_ds,floor([size_h_wdf, size_w_wdf]/wdf_ds));
        
        % udpate patch size
        patch_size_h_new = bottom_right(1) - top_left(1) + 1;
        patch_size_w_new = bottom_right(2) - top_left(2) + 1;
        patch_wdf_size_h_new = bottom_right_wdf(1) - top_left_wdf(1) + 1;
        patch_wdf_size_w_new = bottom_right_wdf(2) - top_left_wdf(2) + 1;
        patch_wdf_size_h_ds = bottom_right_wdf_ds(1) - top_left_wdf_ds(1) + 1;
        patch_wdf_size_w_ds = bottom_right_wdf_ds(2) - top_left_wdf_ds(2) + 1;
        
        
        patch_loc = [top_left; bottom_right];
        patch_wdf_loc = [top_left_wdf; bottom_right_wdf];
        patch_wdf_loc_ds = [top_left_wdf_ds; bottom_right_wdf_ds];
                
        patch_info.location = patch_loc;
        patch_info.size = [patch_size_h_new, patch_size_w_new];
        patch_info.wdf_loc = patch_wdf_loc;
        patch_info.wdf_size = [patch_wdf_size_h_new,patch_wdf_size_w_new];
        patch_info.wdf_loc_ds = patch_wdf_loc_ds;
        patch_info.wdf_ds_size = [patch_wdf_size_h_ds,patch_wdf_size_w_ds];
        patch_info_array_length = length(patch_info_array);
        patch_info_array{patch_info_array_length + 1} = patch_info;
    end
end
end