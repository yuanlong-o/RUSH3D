clc;clear;
close all;

%% this file is the main file for the RUSH3D calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  patched version, for the sake of processing memory and also processing
%  speed.
%  registration module is inserted  after realignment. Dynamic memory calcualtion
%  enabled.


%  last update: 6/29/2021. MW
%  last update: 6/5/2021. YZ
%  last update: 5/23/2021. MW
%  last update: 1/11/2022. MW


%% Path 
% addpath
addpath('background_rejection');
addpath('main demixing');
addpath('preprocess_module');
addpath('reconstruction_module');
addpath('registration');
addpath('seed_generation_module');
addpath(genpath('segmentation_module'));
addpath(genpath('utility'));
%% file path definition

% output dir
main_param.outdir = '..\RUSH3Dresult\1';
% rawdata dir and file name
main_param.rawdata_path = '..\RUSH3Drawdata'; % raw data path (before rotate, resize and rotate);
main_param.first_file_name = 'Rasgrf-Ai148D_1_whiskerstim_3_3x3_50.0ms_Full_Hardware_LaserCount1_211021181554'; % in one video or one capture, the name of first stack
main_param.Nnum = 15;
main_param.Nshift = 3;
main_param.num_rawdata = 800;
main_param.view_range = 5;
main_param.overlap = 21;
main_param.patch_block = [3,4];
main_param.side = [0,15,0,0];
%%
mkdir(main_param.outdir);
rawdata_path = main_param.rawdata_path;
first_file_name = main_param.first_file_name;
input_rawdata_perfix = [rawdata_path, '\\', first_file_name]; % the first file in one capture
view_array = view_config(main_param);
%% Realign parameter
realign_param.num_rawdata = 7600;
realign_param.start_frame = 0; % the number of start frame (the first number is 0, from 0 to N-1)
realign_param.frame_interval = 1; % interval of frame (1 by default: process frame by frame, no jump frame)
realign_param.auto_center_mode = 0; % find the center coordinate of x and y automatically when it is 1 and it will disable center_X and center_Y, otherwise 0
realign_param.conf_name = './utility/realign/3x3.conf.sk.png'; % configuration file for scanning which is corresponding to Nshift
realign_param.center_X = 4000; % the center coordinate x of Light field data (Note that the coordinate is from 0 to N-1 in c++)
realign_param.center_Y = 3001; % the center coordinate y of Light field data (Note that the coordinate is from 0 to N-1 in c++)
realign_param.Nx = 260; % take half number of microlens in x direction (total number: Nx * 2 + 1)
realign_param.Ny = 193; % take half number of microlens in y direction (total number: Nx * 2 + 1)
realign_param.group_mode = 0; % Mode of realign between different frame of rawdata. 0: jump mode (group1: 1-9, group2: 10-18,...); 1: slide window(group1: 1-9, group2: 2-10,...)
realign_param.slight_resize = 0.9991; % slight resize raw data in realign function (1 by default)
realign_param.slight_rotation =  0; % slight rotate raw data in realign function (0 by default) Note that we do not recommend resize and rotate in realign module.
realign_param.rotation =  0; % rotate raw data clockwise (all choice: 0, 90, 180, 270)
realign_param.upsampling_resize = 0; % 1 means resize WDF to 1 * 1, otherwise no resize; 
realign_param.realign_mode = 'LZ'; % realignMode for different scanning sequence (all choices: 'LZ': light path scanning (in RUSH3D). 'ZGX': stage scanning in the opposite direction from 'LZ')
realign_param.bright_scale_normalize = 0; % 1 for normalizing bright scale for each scanning period (for lym) 0 no;
realign_param.skip_zero_frame = 0; % 1 for skip zero frame, 0 no;
realign_param.Nnum = main_param.Nnum; % 15 x 15 pixels behind each microlen. This parameter do not have to change.  
realign_param.Nshift = main_param.Nshift; % scanning for Nshift x Nshift times (all choices: 3: 3 x 3; 5: 5 x 5; 13: 13 x 13)
%% video realign parameter 
video_realign_param = realign_param;
video_realign_param.start_frame = 0;
video_realign_param.upsampling_resize = 1;
video_realign_param.group_mode = 1;
video_realign_param.num_rawdata = main_param.num_rawdata;
%% Realign
%  realign will combine the scanned image 
%  the rotation and resize corrected data
disp('---------------------------Realign---------------------------');

% save folder and file name
outdir = main_param.outdir;
realign_savepath = sprintf('%s\\realign',outdir);
realigndata_name_perfix = sprintf('%s\\realign_%s',realign_savepath,first_file_name);
if ~exist(realign_savepath,'file')
    mkdir(realign_savepath);
end

if realign_param.group_mode == 1
    realign_param.valid_frame_num = realign_param.num_rawdata - realign_param.Nshift^2 - realign_param.start_frame + 1; 
else
    realign_param.valid_frame_num = floor((realign_param.num_rawdata-realign_param.start_frame)/realign_param.Nshift^2);
end

input_rawdata_name = strcat(input_rawdata_perfix,'.0.tiff');
realign_param.auto_center_frame = realign_param.start_frame; % use auto_center_frame to find center automatically under the auto center mode, take the first frame (0 in c++) by default  
realign_param.centerview = strcat(realigndata_name_perfix,'_cv.tiff');

% start realign
tic;
%realign_module(realign_param, input_rawdata_name, realigndata_name_perfix);
t_realign = toc;
fprintf('Realign process done and it takes %.2f secs\n',t_realign);


%% cut patch
% load first wdf

overlap = main_param.overlap;
first_WDF = loadtiff(sprintf('%s_No0.tif', realigndata_name_perfix));
% cut patch
patch_block = main_param.patch_block;
top_left_base = main_param.side(1:2);
top_left = top_left_base * realign_param.Nshift + 1;
right_down_base = [size(first_WDF,1),size(first_WDF,2)]./realign_param.Nshift - main_param.side(3:4);
right_down = [size(first_WDF,1),size(first_WDF,2)] - main_param.side(3:4) * realign_param.Nshift;

% patch preparation
[patch_info_array,wdf_h,wdf_w,wdf_h_ds,wdf_w_ds] = patch_preparation_ovlp_2Nshift(squeeze(first_WDF(:,:,1)), patch_block, top_left_base, right_down_base, ...
                                                             realign_param.Nshift,1,overlap);
%% vessel segmentation
mkdir(sprintf('%s\\vessel_classifier', outdir));
%% Video Realign
%  realign will combine the scanned image 
%  the rotation and resize corrected data
disp('---------------------------video realign---------------------------');

% save folder and file name
video_realign_savepath = sprintf('%s\\video_realign',outdir);
video_realign_name_perfix = sprintf('%s\\video_%s',video_realign_savepath,first_file_name);
if ~exist(video_realign_savepath,'file')
    mkdir(video_realign_savepath);
end

if video_realign_param.group_mode == 1
    video_realign_param.valid_frame_num = video_realign_param.num_rawdata -video_realign_param.Nshift^2 - video_realign_param.start_frame + 1; 
else
    video_realign_param.valid_frame_num = floor((video_realign_param.num_rawdata-video_realign_param.start_frame)/video_realign_param.Nshift^2);
end

input_rawdata_name = strcat(input_rawdata_perfix,'.0.tiff');
video_realign_param.auto_center_frame = video_realign_param.start_frame; % use auto_center_frame to find center automatically under the auto center mode, take the first frame (0 in c++) by default  
video_realign_param.centerview = 'None';

% start realign
tic;
%realign_module(video_realign_param, input_rawdata_name, video_realign_name_perfix);
t_video_realign = toc;
fprintf('video realign process done and it takes %.2f secs\n',t_video_realign);

video_realign_param.Nshift = 1;
video_realign_param.valid_frame_num = video_realign_param.num_rawdata;


%% save parameters
main_param.wdf_h = wdf_h;
main_param.wdf_w = wdf_w;
main_param.wdf_h_ds = wdf_h_ds;
main_param.wdf_h_ds = wdf_w_ds;
main_param.margin_point = [top_left;right_down];
save('.\\param.mat','main_param', '-v7.3');
save(sprintf('%s\\patch_info.mat',outdir),'patch_info_array', '-v7.3');
%% Detrending, Registration, Reconstruction, Calcium parameters
% ------------------- Multiscale detrending parameter ------------------
debg_param.mul_ds_psf = 30; % down sample rate of psf is 12
debg_param.PSF_broader = 10; % cut psf 12 pixel for each side (four sides)
debg_param.psf_end = 101; % psf axial number
debg_param.z_select = 0; % 0 auto select zrange; 1: set manually
debg_param.view_range = 5;
debg_param.writemode = 0; % 0: no output, 2: output last iteration, 1: ouput all result
debg_param.dispmode = 0; % 0: no disp, 1: disp
debg_param.Nbx = 1;
debg_param.Nby = 1; % Block apart 5 x 5 pieces when DAO
debg_param.maxIter = 1 ; % Max iteration times: 2 or 3 is enough
debg_param.AOstar = 1; % 1 for DAO; 0 for no DAO
debg_param.defocus = 1; % 1 for defocus, 0 for no defocus
debg_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
debg_param.margin = 9; % margin overlap
debg_param.first_index = 10;
debg_param.second_index = 18;

debg_param.angle_range = 20;
debg_param.restart = 3;
debg_param.maxframe = 1200;
% ----------------- Registration and std parameter --------------------
reg_param.reg_gSig = 2;
reg_param.reg_bin_width = 200;
reg_param.rankdetrending = 1;
reg_param.maxIter = 10;
reg_param.CoreNum_view = 4;
% ----------------------------- Reconstruction parameter ------------------------------
recon_param.PSF_broader = 50; %% Cut PSF each side for 276 pixels;
recon_param.ds = 3; %% Downsampling rate for PSF
recon_param.psf_end = 101;
recon_param.writemode = 1;
recon_param.dispmode = 1;
recon_param.angle_range = 20;
recon_param.restart = 3;
recon_param.Nbx = 3;
recon_param.Nby = 3; % Block apart 1 x 1 pieces when DAO
recon_param.num_block = 1 ; % Axial block for 10 when forword propagate
recon_param.maxIter = 5 ; % Max iteration times: 2 or 3 is enough
recon_param.view_range = 5; % About 25 Angle views within the iteration
recon_param.AOstar = 1; % 1 for DAO; 0 for no DAO
recon_param.defocus = 1; % 1 for defocus, 0 for no defocus
recon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
recon_param.margin = 9; % margin overlap
recon_param.patchx = 2;
recon_param.patchy = 2;
% --------------------segmentation parameters--------------------
seed_param.pixel_size = 4.6e-6/4.561*recon_param.ds; % this three parameter is not used in reconstruction but use in calcium extraction
seed_param.per_slice_depth = 8e-6; % depth
seed_param.estimate_patch_size = 250;
seed_param.mask3D = 0;
seed_param.neuron_number = 300;
seed_param.neuron_lateral_size = 6;
seed_param.local_constrast_th = 1.45;
seed_param.overlap_ratio = 0.5;
seed_param.boundary_thickness = 2; % for plotwo
seed_param.discard_component_save = false;
seed_param.start_ind = 40;
seed_param.end_ind = 93;
seed_param.max_value = 2000;
seed_param.min_value = 50;
% --------------------iteration related--------------------
demix_param.bg_iter = 10;
demix_param.max_demixing_round = 3;
demix_param.maxIter_NMF = 10;
% --------------------background rejection--------------------
demix_param.oasis_lambda = 0.02;
demix_param.oasis_g = 0.9;
demix_param.lambda_l0 = 0.01;
demix_param.frames_step = 1;
% ------------------- vedio detrending parameter ------------------
viddebg_param = debg_param;
viddebg_param.maxframe = 1000;
% ----------------- Registration and std parameter --------------------
vidreg_param = reg_param;
vidreg_param.reg_bin_width = 80;
vidreg_param.CoreNum_view = 4;
% ----------------- PSF path --------------------
psf_param.correction = 1;
psf_param.debg_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-1000um_to_1000um_zstep_20_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5_ds_30\psf_zmat';
psf_param.recon_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-300um_to_300um_zstep_6_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5_ds_3\psf_zmat';
psf_param.reconpsf_perfix = 'D:\RUSH3Dproject\RUSH3Dpsfblock\20211115';
psf_param.reconpsf_surfix = '_4.561x_lambda_525_axial_-400um_to_400um_zstep_8_NA_0.4_fml_393.3um_OSR_5_ds_3\psf_zmat';
%% save parameters
realign_param.realigndata_name_perfix = realigndata_name_perfix;
video_realign_param.realigndata_name_perfix = video_realign_name_perfix;
save(sprintf('%s\\param_main.mat',outdir),'main_param','realign_param', 'psf_param', 'video_realign_param',...
    'debg_param','reg_param','recon_param','seed_param','demix_param','viddebg_param','vidreg_param','-v7.3');