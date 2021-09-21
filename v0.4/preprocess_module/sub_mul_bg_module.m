function curr_video = sub_mul_bg_module(psf,debgrecon_param,curr_video,vid_recon_mul_savepath,vid_debg_mul_savepath,filename,frame,channel)
%% This program is used for substrcating multiscale

Nnum = debgrecon_param.Nnum;
Nshift = debgrecon_param.Nshift;
mul_ds_psf = debgrecon_param.mul_ds_psf;
view_array = debgrecon_param.view_array;


curr_wdf = single(curr_video);
%% reconstruct with normal deconvolution

gpuDevice;
vid_recon_mul_savename = strcat('vid_', filename); % share the same file name
vid_recon_name_perfix = strcat(vid_recon_mul_savepath, '\\', vid_recon_mul_savename);

upsampling = Nnum / Nshift / mul_ds_psf;

vid_WDF_up = imresize(curr_wdf, ...
    [floor(size(curr_wdf,1)*upsampling/2)*2+1,floor(size(curr_wdf,2)*upsampling/2)*2+1]);

vid_volume = ones(size(vid_WDF_up,1),size(vid_WDF_up,2),size(psf,4));
vid_volume = vid_volume./size(vid_volume,3);
vid_volume = reconstruction_module(debgrecon_param, vid_volume, psf, vid_WDF_up, vid_recon_name_perfix, view_array, frame);

%% substract background

vid_volume_z_side = vid_volume;
vid_volume_z_side(:,:,debgrecon_param.first_index: debgrecon_param.second_index) = 0;
vid_volume_z_mid = vid_volume;
vid_volume_z_mid(:,:,[1:debgrecon_param.first_index-1, debgrecon_param.second_index+1:end]) = 0;

wdf_side = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),length(view_array),'single');

global bg_ratio

wdf_debg = zeros(size(curr_wdf),'single');

if frame-channel == 0
    wdf_mid = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),length(view_array),'single');
    for v = 1: length(view_array)
        wdf_side(:,:,v) = forwardProjectSingleView(squeeze(psf(:,:,v,:)),vid_volume_z_side);
        wdf_mid(:,:,v)  = forwardProjectSingleView(squeeze(psf(:,:,v,:)),vid_volume_z_mid);
        
        bg_ratio(v,channel+1) = mean(mean(wdf_side(:,:,v)))/mean(mean(wdf_mid(:,:,v)));
        bg_ratio(v,channel+1) =mean(mean(curr_wdf(:,:,v)))/mean(mean(wdf_side(:,:,v)))*bg_ratio(v)/(1+bg_ratio(v));
    end
end

for v = 1: length(view_array)
    wdf_side(:,:,v) = forwardProjectSingleView(squeeze(psf(:,:,v,:)),vid_volume_z_side);
    tmp = wdf_side(:,:,v) * bg_ratio(v,channel+1);
    tmp = gather(tmp);
    tmp = imresize(tmp,[size(curr_wdf,1),size(curr_wdf,2)]);
    wdf_debg(:,:,v) = curr_wdf(:,:,v) - tmp;
end

wdf_debg = max(wdf_debg,0);

curr_video = wdf_debg;

% saveastiff(uint16(curr_video / max(curr_video(:))*65535), ...
%     sprintf('%s\\wdf_debg_vid_c%d_%d.tif', vid_debg_mul_savepath,channel,frame-channel));

%%
end