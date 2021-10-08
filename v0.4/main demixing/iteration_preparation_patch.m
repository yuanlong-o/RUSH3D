function processed_video = iteration_preparation_patch( savegroup, reg_path, valid_frame, ...
                                                        specify_wigner, curr_patch_info, ...
                                                        psf_param, seed_param, S_init)
%  load wigners, and subtract marginal view information from central view
%  output is the processed video.
%  last update: 6/5/2021. YZ

%% take into account the scanning issue

vol_ds = seed_param.vol_ds;
wdf_ds = seed_param.wdf_ds;
NN_ds = seed_param.NN_ds;
ds_ratio = NN_ds * wdf_ds / vol_ds;
side = psf_param.side;


buf_S = S_init{1, 1};
buf_S = full(buf_S);
[size_h, size_w] = size(buf_S);

% query the overall size


% buf_wigner_video = loadtiff(sprintf('%s\\reg_view_%d_g%d.tiff', reg_path, specify_wigner(1),1));
% buf_wigner_video_cut = buf_wigner_video(side+1:end-side,side+1:end-side,:);
% buf_wigner_video_ds = imresize(buf_wigner_video_cut,1/ds_ratio);

% [size_h_global, size_w_global, ~] = size(buf_wigner_video_ds);

% % determine cut position for h
% cut_start_h = max(ceil(curr_patch_info.location(1, 1) / ds_ratio), 1);
% if cut_start_h + size_h - 1 > size_h_global % count the size from the end
%     cut_end_h = size_h_global;
%     cut_start_h = cut_end_h - size_h + 1;
% else
%     cut_end_h  = cut_start_h + size_h - 1;
% end
% 
% % the same for w
% cut_start_w = max(ceil(curr_patch_info.location(1, 2) / ds_ratio), 1);
% if cut_start_w + size_w - 1 > size_w_global % count the size from the end
%     cut_end_w = size_w_global;
%     cut_start_w = cut_end_w - size_w + 1;
% else
%     cut_end_w  = cut_start_w + size_w - 1;
% end

% make sure the size is the same as S
num_wigner = size(specify_wigner, 1);
processed_video = zeros(size_h, size_w, valid_frame, num_wigner, 'single');
%% main load module
for j = 1 : num_wigner % number of specified wigner
    curr_video = [];
    for g_id = 1 : savegroup
        curr_video = cat(3,curr_video,loadtiff(sprintf('%s\\reg_view_%d_g%d.tiff', reg_path, specify_wigner(j),g_id)));
    end
    curr_video_cut = curr_video(side+1:end-side,side+1:end-side,:);
    curr_video_cut = single(curr_video_cut) / 65535;
    curr_video_ds = imresize(curr_video_cut,1/wdf_ds);
    processed_video(:, :, :, j) = curr_video_ds(curr_patch_info.wdf_loc(1, 1) : curr_patch_info.wdf_loc(2, 1), ...
        curr_patch_info.wdf_loc(1, 2) : curr_patch_info.wdf_loc(2, 2),:);
end 



end