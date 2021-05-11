function Xguess = reconstruction_module(param, wdfname, psfpath)
%% this function run reconstructions for RUSH3D
%  wdfname: summarized raw signals
%  psfpath: path that contains psfs
%  last update: 5/11/2021, YZ

%% parameters querying
gpuDevice_id = param.gpuDevice_id;
M = param.M; %% Magnification = 3.17
Nshift = param.Nshift; %% Scanning Times = 3 x 3
Nnum = param.Nnum; %% 15 x 15 pixels behind each MicroLen
PSF_broader = param.PSF_broader; %% Cut PSF each side for 276 pixels;
PSF_size = param.PSF_size; %%��x,y,u,v,z��---->(x,y,u,z) no meaning
Block_num = param.Block_num; %% no meaning

% ------------------------- Reconstruction parameter ------------------------- % 

Angle_range = param.Angle_range; %% About 25 Angle views within the iteration
AOstar = param.AOstar; %% 1 for DAO; 0 for no DAO 
maxIter = param.maxIter; %% Max iteration times: 2 or 3 is enough
rotWDF = param.rotWDF; %% no meaning
DownSampling_for_PSF = param.DownSampling_for_PSF; %% Downsampling rate for PSF
UpSamplingRate_for_WDF = Nnum/Nshift/DownSampling_for_PSF; %% UpsamlpingRate for WDF, which is corresponding to the PSF Size.
% psf_layer = [1:2:51,52:1:151,153:2:201]; % PSF to be loaded []

% ------------------------- DAO parameter ---------------------------

defocus = param.defocus;  %% 1 for defocus, 0 for no defocus
Nbx = param.Nbx;
Nby = param.Nby; %% Block apart 5 x 5 pieces when DAO 
num_block = param.num_block;    %% Axial block for 10 when forword propagate 
threshhold = param.threshhold;  %% Shift should not be allowed to exceed [-25,25]
sidelobe = param.sidelobe; %% Each 

% Load WDF
tmp = double(imread(wdfname,1));
WDF = zeros(size(tmp,1),size(tmp,2),Nnum,Nnum);


gpuDevice(gpuDevice_id)

for u = 1:Nnum
    for v = 1: Nnum
        k = (u-1)*Nnum + v;
        tmp = double(imread(wdfname,k));
        tmp(tmp<0)=0;
        WDF(:,:,u,v) = tmp;
        disp(['u = ', num2str(u), ', v = ', num2str(v), 'has been load !']);
    end
end

%% load PSF
load([psfpath,'/psf_sim_',num2str(M),'_',num2str(1),'.mat']);
psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:,:);
psf_z = imresize(psf_z,[floor(size(psf_z,1) / DownSampling_for_PSF),floor(size(psf_z,2) / DownSampling_for_PSF)],'cubic');
psf = zeros(size(psf_z,1),size(psf_z,2),size(psf_z,3),size(psf_z,4),size(psf_layer,2));
layer_id = 1;
for layer_num = psf_layer
    load([psfpath,'\psf_sim_',num2str(M),'_',num2str(layer_num),'.mat']);
    psf_z = psf_z(PSF_broader:end-PSF_broader+1,PSF_broader:end-PSF_broader+1,:,:,:);
    psf_z = imresize(psf_z,[floor(size(psf_z,1) / DownSampling_for_PSF),floor(size(psf_z,2) / DownSampling_for_PSF)],'cubic');
    if layer_num < 51 || layer_num > 151
        psf_z = 2*psf_z;
    elseif layer_num == 51 || layer_num == 151
        psf_z = 1.5*psf_z;
    end
    psf(:,:,:,:,layer_id) = psf_z;
    disp(['layer_id = ', num2str(layer_id), '  has been load !']);
    layer_id = layer_id + 1;
end


%% Reconstruction
WDF_resize = imresize(WDF,[floor(size(WDF,1) * UpSamplingRate_for_WDF),floor(size(WDF,2) * UpSamplingRate_for_WDF)],'cubic');
WDF_resize = rot90(WDF_resize, 2*rotWDF);

Xguess=ones(size(WDF_resize,1),size(WDF_resize,2),size(psf,5));
Xguess=Xguess./sum(Xguess(:)).*sum(Xguess(:))./(size(Xguess,3)*size(Xguess,4));
Xguess = deconvRLGPUwithAO_wmr(maxIter,Xguess,WDF_resize,psf,Nnum,0, [],AOstar,Angle_range,Block_num, defocus, Nbx,Nby,num_block,threshhold,sidelobe);
%imwriteTFSK(single(Xguess),[recondata_folder3,num2str(Nshift),'_angle=',num2str(Angle_range),'_AO=',num2str(AOstar),'_normal','.png']);


end