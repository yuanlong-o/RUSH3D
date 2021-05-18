function [psf,psf_param] = psf_load_module(psf_param, psfpath)
%% PSF Load Module Load PSF
% This program is used to load psf with multiscale and downsampling.
% Last update: 05/15/2021. MW 
% Last update: 05/18/2021. YZ


% Input:
% psf_param
% psf_path       psf_path_perfix

% Output:
% psf_param.psfsize
% psf

M = psf_param.M; % Magnification = 3.17
PSF_broader = psf_param.PSF_broader; % Cut PSF each side for 276 pixels;
downsampling_rate = psf_param.downsampling_rate; % Downsampling rate for PSF
psf_layer_position = psf_param.psf_layer_position; % multiscale layer load

layer = psf_layer_position(1);
psf_z = load([psfpath ,'/psf_sim_',num2str(M),'_',num2str(layer),'.mat']);
if isstruct(psf_z )
   psf_z  = psf_z.psf_z; 
end

psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:,:);
psf_z = imresize(psf_z,downsampling_rate);
psf = zeros(size(psf_z,1),size(psf_z,2),size(psf_z,3),size(psf_z,4),size(psf_layer_position,2));
psfsize = size(psf);

count = 1;
for layer = psf_layer_position
    psf_z = load([psfpath ,'/psf_sim_',num2str(M),'_',num2str(layer),'.mat']);
    % handle the structure
    if isstruct(psf_z )
       psf_z  = psf_z.psf_z; 
    end
    psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:,:);
    psf_z = imresize(psf_z,downsampling_rate);
    if layer < 51 || layer > 151
        psf_z = 2*psf_z;
    elseif layer == 51 || layer == 151
        psf_z = 1.5*psf_z;
    end
    psf(:,:,:,:,count) = psf_z;
    count = count+1;
end
psf_param.psfsize = psfsize;
end