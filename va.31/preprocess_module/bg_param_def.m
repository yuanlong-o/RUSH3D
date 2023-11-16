function debgrecon_param = bg_param_def(psf,debgrecon_param,curr_video,vid_recon_mul_savepath,filename)

%% This program is used for substrcating multiscale

Nnum = debgrecon_param.Nnum;
Nshift = debgrecon_param.Nshift;
mul_ds_psf = debgrecon_param.mul_ds_psf;
view_array = debgrecon_param.view_array;
z_select = debgrecon_param.z_select;
temp_writemode = debgrecon_param.writemode;
debgrecon_param.writemode = 2;
first_index = debgrecon_param.first_index;
second_index = debgrecon_param.second_index;


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
vid_volume = recon_module(debgrecon_param, vid_volume, psf, vid_WDF_up, vid_recon_name_perfix, view_array, 0);
%% calculate z range
if z_select == 0
    z_energy = squeeze(mean(vid_volume,[1,2]));
    [max_pk_l,zloc] = findpeaks(z_energy);
    min_pk = min(z_energy);
    if ~isempty(zloc)
        max_pk = squeeze(mean(max_pk_l));
    else
        max_pk = max(z_energy);
    end
    mid_th = (min_pk + max_pk)/2;
    z_list = find(z_energy>mid_th);
    z_list((z_list - second_index) > (size(psf,4)/5)) = []; 
    z_list((first_index - z_list) > (size(psf,4)/5)) = []; 
    debgrecon_param.first_index = min(z_list)-1;
    debgrecon_param.second_index = max(z_list)+1;
end
%% calculate bg_ratio

bg_ratio = gpuArray.zeros(length(view_array),1);
vid_volume_z_side = vid_volume;
vid_volume_z_side(:,:,debgrecon_param.first_index: debgrecon_param.second_index) = 0;
vid_volume_z_mid = vid_volume;
vid_volume_z_mid(:,:,[1:debgrecon_param.first_index-1, debgrecon_param.second_index+1:end]) = 0;

wdf_side = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),length(view_array),'single');
wdf_mid = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),length(view_array),'single');
for v = 1: length(view_array)
    wdf_side(:,:,v) = forwardProjectSingleView(squeeze(psf(:,:,v,:)),vid_volume_z_side);
    wdf_mid(:,:,v)  = forwardProjectSingleView(squeeze(psf(:,:,v,:)),vid_volume_z_mid);
    
    bg_ratio(v) = mean(mean(wdf_side(:,:,v)))/mean(mean(wdf_mid(:,:,v)));
    bg_ratio(v) =mean(mean(curr_wdf(:,:,v)))/mean(mean(wdf_side(:,:,v)))*bg_ratio(v)/(1+bg_ratio(v));
end

bg_ratio = gather(bg_ratio);

debgrecon_param.bg_ratio = bg_ratio;
debgrecon_param.writemode = temp_writemode;
% saveastiff(uint16(curr_video / max(curr_video(:))*65535), ...
%     sprintf('%s\\wdf_debg_vid_c%d_%d.tif', vid_debg_mul_savepath,channel,frame-channel));
%%
end

 