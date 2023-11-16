function curr_video = sub_mul_bg_module(psf,debgrecon_param,curr_video,vid_recon_mul_savepath,vid_debg_mul_savepath,filename,frame,channel)
%% This program is used for substrcating multiscale

Nnum = debgrecon_param.Nnum;
Nshift = debgrecon_param.Nshift;
mul_ds_psf = debgrecon_param.mul_ds_psf;

curr_wdf = reshape(curr_video,size(curr_video,1),size(curr_video,2),Nnum,Nnum);
curr_wdf = permute(curr_wdf,[1,2,4,3]);

%% reconstruct with normal deconvolution

gpuDevice;
vid_recon_mul_savename = strcat('vid_', filename); % share the same file name
vid_recon_name_perfix = strcat(vid_recon_mul_savepath, '\\', vid_recon_mul_savename);

upsampling = Nnum / Nshift / mul_ds_psf;

vid_WDF_up = imresize(curr_wdf, ...
    [floor(size(curr_wdf,1)*upsampling/2)*2+1,floor(size(curr_wdf,2)*upsampling/2)*2+1],'cubic');
vid_volume = ones(size(vid_WDF_up,1),size(vid_WDF_up,2),size(psf,5));
vid_volume = vid_volume./sum(vid_volume(:)).*sum(vid_volume(:))./(size(vid_volume,3)*size(vid_volume,4));
vid_volume = reconstruction_module(debgrecon_param, vid_volume, psf, vid_WDF_up, vid_recon_name_perfix, frame);

%% substract background

vid_volume_z_side = vid_volume;
vid_volume_z_side(:,:,debgrecon_param.first_index+1: debgrecon_param.second_index) = 0;
vid_volume_z_mid = vid_volume;
vid_volume_z_mid(:,:,[1:debgrecon_param.first_index, debgrecon_param.second_index+1:end]) = 0;

wdf_side = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),debgrecon_param.Nnum,debgrecon_param.Nnum,'single');

global bg_ratio

wdf_debg = gpuArray.zeros(size(curr_wdf),'single');

if frame-channel == 0
    wdf_mid = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),debgrecon_param.Nnum,debgrecon_param.Nnum,'single');
    for u = 1:debgrecon_param.Nnum
        for v = 1:debgrecon_param.Nnum
            wdf_side(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_side);
            wdf_mid(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_mid);
        
            bg_ratio(u,v,channel+1) = mean(mean(wdf_side(:,:,u,v)))/mean(mean(wdf_mid(:,:,u,v)));
            bg_ratio(u,v,channel+1) =mean(mean(curr_wdf(:,:,u,v)))/mean(mean(wdf_side(:,:,u,v)))*bg_ratio(u,v)/(1+bg_ratio(u,v));

        end
    end
end

    
for u = 1:debgrecon_param.Nnum
    for v = 1:debgrecon_param.Nnum
        wdf_side(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_side);

        tmp = wdf_side(:,:,u,v) * bg_ratio(u,v,channel+1);
        tmp = imresize(tmp,[size(curr_wdf,1),size(curr_wdf,2)]);

        wdf_debg(:,:,u,v) = curr_wdf(:,:,u,v) - tmp;
    end
end
wdf_debg = max(wdf_debg,0);

wdf_debg = permute(wdf_debg,[1,2,4,3]);
wdf_debg = gather(wdf_debg);
curr_video = reshape(wdf_debg,size(wdf_debg,1),size(wdf_debg,2),Nnum*Nnum);
saveastiff(uint16(curr_video / max(curr_video(:))*65535), ...
    sprintf('%s\\wdf_debg_vid_c%d_%d.tif', vid_debg_mul_savepath,channel,frame-channel));

% Output=uint16(wdf_side./max(wdf_side(:)).*65535);
% for idu=1:psf_param.Nnum
%     for idv = 1:psf_param.Nnum
%         ii = (idu-1)*psf_param.Nnum+idv;
%         if ii==1
%             imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(vid_debg_mul_savepath,'\\','wdf_side_vid',num2str(frame),'.tif'));
%         else
%             imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(vid_debg_mul_savepath,'\\','wdf_side_vid',num2str(frame),'.tif'),'WriteMode', 'append');
%         end
%     end
% end

% Output=uint16(wdf_debg./max(wdf_debg(:)).*65535);
% for idu=1:psf_param.Nnum
%     for idv = 1:psf_param.Nnum
%         ii = (idu-1)*psf_param.Nnum+idv;
%         if ii==1
%             imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(vid_debg_mul_savepath,'\\','wdf_debg_vid',num2str(frame),'.tif'));
%         else
%             imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(vid_debg_mul_savepath,'\\','wdf_debg_vid',num2str(frame),'.tif'),'WriteMode', 'append');
%         end
%     end
% end
%%
end