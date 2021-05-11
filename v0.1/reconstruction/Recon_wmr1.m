clc
clear
close all
gpuDevice(1);

%% Parameters

% ------------------------- Main parameter --------------------------- %

M = 3.1746; %% Magnification = 3.17
Nshift = 3; %% Scanning Times = 3 x 3
Nnum = 15; %% 15 x 15 pixels behind each MicroLen
PSF_broader = 276; %% Cut PSF each side for 276 pixels;
PSF_size = 5; %%（x,y,u,v,z）---->(x,y,u,z) no meaning
Block_num = 0; %% no meaning


% ------------------------- Reconstruction parameter ------------------------- % 

Angle_range = 25; %% About 25 Angle views within the iteration
AOstar = 1; %% 1 for DAO; 0 for no DAO 
maxIter = 3; %% Max iteration times: 2 or 3 is enough
rotWDF = 0; %% no meaning
DownSampling_for_PSF = 1; %% Downsampling rate for PSF
UpSamplingRate_for_WDF = Nnum/Nshift/DownSampling_for_PSF; %% UpsamlpingRate for WDF, which is corresponding to the PSF Size.
psf_layer = [1:2:51,52:1:151,153:2:201]; % PSF to be loaded []

% ------------------------- DAO parameter ---------------------------

defocus = 1; %% 1 for defocus, 0 for no defocus
Nbx = 5; Nby = 5; %% Block apart 5 x 5 pieces when DAO 
num_block= 18 ;    %% Axial block for 10 when forword propagate 
threshhold = 25;  %% Shift should not be allowed to exceed [-25,25]
sidelobe=9; %% Each 

%% Wigner Path

rawdata_folder1 = 'Z:/数据采集/';
rawdata_folder2 = 'gyd0309/';
rawdata_folder12 = strcat(rawdata_folder1, rawdata_folder2);
rawdatafilename_perfix = 'X1_NORMAL_3x3_25.0ms_HBar_Hardware_LaserCount1_210309202631';

wdfdata_folder1 = 'D:/RichardW2021/RUSH3DResult/Preprocess/';
wdfdata_folder2 = rawdata_folder2;
wdfdata_folder12 = strcat(wdfdata_folder1, wdfdata_folder2);
wdfdata_folder3 = strcat(wdfdata_folder12, rawdatafilename_perfix,'/');
wdfdata_folder4 = strcat(wdfdata_folder3, 'WDF/');

wdfname = [wdfdata_folder4,'rank1__No0.tif'];

%% PSF Path

psfpath = 'Z:/PSF/psf_zww/20200125_genepsf_3.1746x_sim_neg400T400_dz4_15Nnum_OSR5';

%% Save Path

recondata_folder1= 'D:/RichardW2021/RUSH3DResult/Recon/';
recondata_folder2 = rawdata_folder2;
recondata_folder12 = strcat(recondata_folder1, recondata_folder2);
recondata_folder3 = [recondata_folder12,rawdatafilename_perfix,'/'];

if ~exist(recondata_folder3,'file')
    mkdir(recondata_folder3);
end

%% Load WDF

tmp = double(imread(wdfname,1));
WDF = zeros(size(tmp,1),size(tmp,2),Nnum,Nnum);

for u = 1:Nnum
    for v = 1: Nnum
        k = (u-1)*Nnum + v;
        tmp = double(imread(wdfname,k));
        tmp(tmp<0)=0;
        WDF(:,:,u,v) = tmp;
        disp(['u = ', num2str(u), ', v = ', num2str(v), 'has been load !']);
    end
end

%%  Load PSF


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
 Xguess = deconvRLGPUwithAO_wmr(maxIter,Xguess,WDF_resize,psf,Nnum,0,recondata_folder3,AOstar,Angle_range,Block_num, defocus, Nbx,Nby,num_block,threshhold,sidelobe);
%imwriteTFSK(single(Xguess),[recondata_folder3,num2str(Nshift),'_angle=',num2str(Angle_range),'_AO=',num2str(AOstar),'_normal','.png']);
