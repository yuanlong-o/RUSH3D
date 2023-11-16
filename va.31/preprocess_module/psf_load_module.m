function [psf,psf_param] = psf_load_module(psf_param, psfpath, view_array)
%% PSF Load Module Load PSF
% This program is used to load psf with multiscale and downsampling.
% Last update: 05/15/2021. MW 
% Last update: 05/18/2021. YZ
% Last update: 09/16/2021. YZ

% Input:
% psf_param
% psf_path       psf_path_perfix

% Output:
% psf_param.psfsize
% psf
PSF_broader = psf_param.PSF_broader; 
psf_layer_position = psf_param.psf_layer_position; % multiscale layer load
layer = psf_layer_position(1);
load(sprintf('%s\\psf_sim%d.mat',psfpath,layer),'psf_z'),

psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:);
psf = zeros(size(psf_z,1),size(psf_z,2),length(view_array),size(psf_layer_position,2));
psfsize = size(psf);

count = 1;
for layer = psf_layer_position
    load(sprintf('%s\\psf_sim%d.mat',psfpath,layer),'psf_z'),
    psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:);
    for v = 1 : length(view_array)
        psf(:,:,v,count) = psf_z(:,:,view_array{v}(1),view_array{v}(2));
    end
    count = count + 1;
end
psf_param.psfsize = psfsize;
end