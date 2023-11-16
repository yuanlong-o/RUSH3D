function processed_video = video_preparation_patch( savegroup, reg_path, valid_frame, ...
                                                        specify_wigner, curr_patch_info, ...
                                                        psf_param, seed_param, S_init,shift_block)
%  load wigners, and subtract marginal view information from central view
%  output is the processed video.
%  last update: 6/5/2021. YZ

%% take into account the scanning issue

buf_S = S_init{1, 1};
buf_S = full(buf_S);
[size_h, size_w] = size(buf_S);

num_wigner = size(specify_wigner, 1);
processed_video = zeros(size_h, size_w, valid_frame, num_wigner, 'single');
%% main load module

for j = 1 : num_wigner % number of specified wigner
    curr_video = [];

    for g_id = 1 : savegroup
        curr_video = cat(3,curr_video,loadtiffsub(sprintf('%s\\vidreg_view_%d_g_%d.tiff', reg_path, specify_wigner(j),g_id),...
            curr_patch_info.wdf_loc_ds(1,:)+shift_block(:,j)',curr_patch_info.wdf_loc_ds(2,:)+shift_block(:,j)'));
    end
    curr_video = single(curr_video) / 65535;
    processed_video(:, :, :, j) = curr_video;
    clear curr_video;
end 



end