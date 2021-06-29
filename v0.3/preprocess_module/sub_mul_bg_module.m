function curr_video = sub_mul_bg_module(psf,psf_param,preprocess_param,curr_video,vid_recon_mul_savepath,vid_debg_mul_savepath,filename,frame)

%% This program is used for substrcating multiscale

Nnum = preprocess_param.Nnum;
mul_ds_psf = preprocess_param.mul_ds_psf;
Nshift = preprocess_param.Nshift;
rotWDF = preprocess_param.rotWDF;

curr_wdf = reshape(curr_video,size(curr_video,1),size(curr_video,2),Nnum,Nnum);
curr_wdf = permute(curr_wdf,[1,2,4,3]);
%% reconstruction param

recon_param.writemode = 0; % 0: no output, 2: output last iteration, 1: ouput all result
recon_param.dispmode = 0; % 0: no disp, 1: disp
recon_param.Nbx = 1;
recon_param.Nby = 1; % Block apart 5 x 5 pieces when DAO

recon_param.maxIter = 1 ; % Max iteration times: 2 or 3 is enough
recon_param.angle_range = 12; % About 25 Angle views within the iteration

recon_param.AOstar = 1; % 1 for DAO; 0 for no DAO
recon_param.defocus = 1; % 1 for defocus, 0 for no defocus
recon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
recon_param.margin = 9; % margin overlap
recon_param.estimate_patch_size = 700;

%% reconstruct with normal deconvolution

gpuDevice;

vid_recon_mul_savename = strcat('vid_', filename); % share the same file name
vid_recon_name_perfix = strcat(vid_recon_mul_savepath, '\\', vid_recon_mul_savename);

upsampling = Nnum / Nshift / mul_ds_psf;

vid_WDF_up = rot90(imresize(curr_wdf, ...
    [floor(size(curr_wdf,1)*upsampling/2)*2+1,floor(size(curr_wdf,2)*upsampling/2)*2+1],'cubic'),2*rotWDF);
vid_volume = ones(size(vid_WDF_up,1),size(vid_WDF_up,2),size(psf,5));
vid_volume = vid_volume./sum(vid_volume(:)).*sum(vid_volume(:))./(size(vid_volume,3)*size(vid_volume,4));
vid_volume = reconstruction_module(psf_param, recon_param, vid_volume, psf, vid_WDF_up, vid_recon_name_perfix, frame);
%%

vid_volume_z_side = vid_volume;
vid_volume_z_side(:,:,psf_param.first_index+1: psf_param.second_index) = 0;

vid_volume_z_mid = vid_volume;
vid_volume_z_mid(:,:,[1:psf_param.first_index, psf_param.second_index+1:end]) = 0;

wdf_side = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),psf_param.Nnum,psf_param.Nnum,'single');
% wdf_mid = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),psf_param.Nnum,psf_param.Nnum,'single');
global bg_ratio
% bg_ratio = gpuArray.zeros(psf_param.Nnum,psf_param.Nnum,'single');
% shift_x = gpuArray.zeros(psf_param.Nnum,psf_param.Nnum,'single');
% shift_y = gpuArray.zeros(psf_param.Nnum,psf_param.Nnum,'single');

wdf_debg = gpuArray.zeros(size(curr_wdf),'single');
% wdf_mid_u = gpuArray.zeros(size(curr_wdf),'single');
% [coor_x,coor_y] = meshgrid(1:size(curr_wdf,2),1:size(curr_wdf,1));
if frame == 0
    wdf_mid = gpuArray.zeros(size(vid_volume,1),size(vid_volume,2),psf_param.Nnum,psf_param.Nnum,'single');
    for u = 1:psf_param.Nnum
        for v = 1:psf_param.Nnum
            wdf_side(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_side);
            wdf_mid(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_mid);
        
            bg_ratio(u,v) = mean(mean(wdf_side(:,:,u,v)))/mean(mean(wdf_mid(:,:,u,v)));
            bg_ratio(u,v) =mean(mean(curr_wdf(:,:,u,v)))/mean(mean(wdf_side(:,:,u,v)))*bg_ratio(u,v)/(1+bg_ratio(u,v));

        end
    end
end

    
for u = 1:psf_param.Nnum
    for v = 1:psf_param.Nnum
        wdf_side(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_side);
%         wdf_mid(:,:,u,v) = forwardProjectSingleView(squeeze(psf(:,:,u,v,:)),vid_volume_z_mid);

%         wdf_mid_u(:,:,u,v) = imresize(wdf_mid(:,:,u,v),[size(curr_wdf,1),size(curr_wdf,2)]);
%         
%         corr_map=gather(normxcorr2(squeeze(curr_wdf(:,:,u,v)),squeeze(wdf_mid_u(:,:,u,v))));
%         [shift_x(u,v),shift_y(u,v)] = find(corr_map == max(corr_map(:)));
%         shift_x(u,v) = shift_x(u,v) - size(curr_wdf,1);
%         shift_y(u,v) = shift_y(u,v) - size(curr_wdf,2);
%         if abs(shift_x(u,v)) > 10
%             shift_x(u,v) = 0;
%         end
%         if abs(shift_y(u,v)) > 10
%             shift_y(u,v) = 0;
%         end
%         bg_ratio(u,v) = mean(mean(wdf_side(:,:,u,v)))/mean(mean(wdf_mid(:,:,u,v)));
        tmp = wdf_side(:,:,u,v) * bg_ratio(u,v);
        tmp = imresize(tmp,[size(curr_wdf,1),size(curr_wdf,2)]);
%         tmp = interp2(coor_x,coor_y,tmp,coor_x+shift_y(u,v),coor_y+shift_x(u,v));
        wdf_debg(:,:,u,v) = curr_wdf(:,:,u,v) - tmp;
    end
end
wdf_debg = max(wdf_debg,0);

%wdf_side = gather(wdf_side);
%wdf_mid = gather(wdf_mid);
%wdf_mid_u = gather(wdf_mid_u);
wdf_debg = permute(wdf_debg,[1,2,4,3]);
wdf_debg = gather(wdf_debg);
curr_video = reshape(wdf_debg,size(wdf_debg,1),size(wdf_debg,2),Nnum*Nnum);
saveastiff(uint16(curr_video / max(curr_video(:))*65535), ...
    sprintf('%s\\wdf_debg_vid_%d.tif', vid_debg_mul_savepath,frame));

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