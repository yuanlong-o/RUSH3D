function processed_video = iteration_preparation_patch( savegroup, reg_path, valid_frame, ...
                                                        specify_wigner, curr_patch_info, ...
                                                        psf_param, seed_param, S_init)
%  load wigners, and subtract marginal view information from central view
%  output is the processed video.
%  last update: 6/5/2021. YZ

%% take into account the scanning issue


outdir = seed_param.outdir;


buf_S = S_init{1, 1};
buf_S = full(buf_S);
[size_h, size_w] = size(buf_S);

num_wigner = size(specify_wigner, 1);
processed_video = zeros(size_h, size_w, valid_frame, num_wigner, 'single');
%% main load module
if psf_param.shiftflag == 1
    load(sprintf('%s\\shiftmapblock.mat',outdir),'map_wavshape');
end
for j = 1 : num_wigner % number of specified wigner
    curr_video = [];

    for g_id = 1 : savegroup
        curr_video = cat(3,curr_video,loadtiffsub(sprintf('%s\\reg_view_%d_g_%d.tiff', reg_path, specify_wigner(j),g_id),...
            curr_patch_info.wdf_loc(1,:),curr_patch_info.wdf_loc(2,:)));
    end
    curr_video = single(curr_video) / 65535;
    if psf_param.shiftflag == 1
        map_wavshape_cut = squeeze(map_wavshape(curr_patch_info.wdf_loc(1,1):curr_patch_info.wdf_loc(2,1),...
            curr_patch_info.wdf_loc(1,2):curr_patch_info.wdf_loc(2,2),specify_wigner(j),:));
        map_wavshape_cut = imresize(map_wavshape_cut,[curr_patch_info.wdf_ds_size(1),curr_patch_info.wdf_ds_size(2)]);
        [coordinate1,coordinate2]=meshgrid(1:size(map_wavshape_cut,2),1:size(map_wavshape_cut,1));
        curr_video = imresize(curr_video,[curr_patch_info.wdf_ds_size(1),curr_patch_info.wdf_ds_size(2)]);
        for f = 1: size(curr_video,3)
            processed_video(:, :, f, j) = interp2(coordinate1,coordinate2,squeeze(curr_video(:,:,f)),...
                coordinate1+map_wavshape_cut(:,:,2),coordinate2+map_wavshape_cut(:,:,1),'cubic',0);
        end
    else
        processed_video(:, :, :, j) = imresize(curr_video,[curr_patch_info.wdf_ds_size(1),curr_patch_info.wdf_ds_size(2)]);
    end

    clear curr_video;
end 



end