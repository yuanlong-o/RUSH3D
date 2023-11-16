clc;clear;
close all;

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  patched version, for the sake of processing memory and also processing
%  speed.
%  registration module is inserted  after realignment. Dynamic memory calcualtion
%  enabled.

%% addpath
addpath('../background_rejection');
addpath('../main demixing');
addpath('../preprocess_module');
addpath('../reconstruction_module');
addpath('../registration');
addpath('../seed_generation_module');
addpath(genpath('../segmentation_module'));
addpath(genpath('../utility'));

%% colormap

color_scheme_npg = [...
    0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
         0    0.6275    0.5294; ...
    0.2353    0.3294    0.5333; ...
    0.9529    0.6078    0.4980; ...
    0.5176    0.5686    0.7059; ...
    0.5686    0.8196    0.7608; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];

%% file path definition

% output dir
% outdir = 'D:\RUSH3Dproject\RUSH3Dresult\0816\5x_br_tb_0.00_bg_0.00_lfm\test_patch'; % MW path
outdir = 'D:\RUSH3Dproject\RUSH3Dresult\20210928\rasgrf-ai148d_5xTMP_whisker3\test'; % YZ
mkdir(outdir)

% rawdata dir and file name
% rawdata_path = 'D:\RUSH3Dproject\RUSH3Drawdata\0816\5x_br_tb_0.00_bg_0.00_lfm_uni_N_975_stack'; % raw data path (before rotate, resize and rotate);
rawdata_path = 'D:\RUSH3Dproject\RUSH3Drawdata\20210928\rasgrf-ai148d_5xTMP_whisker3'; % raw data path (before rotate, resize and rotate);

first_file_name = 'rasgrf-ai148d_5xTMP_whisker3_3x3_50.0ms_Full_Hardware_LaserCount1_210928213137'; % in one video or one capture, the name of first stack
input_rawdata_perfix = [rawdata_path, '\\', first_file_name]; % the first file in one capture

% debg psfpath
% debg_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\10x_lambda_488_axial_-600um_to_600um_zstep_20_NA_0.5_fml_875um_ftl_265mm_OSR_5_ds_12\psf_zmat';
debg_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-1000um_to_1000um_zstep_20_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5_ds_30\psf_zmat';
recon_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\09265.761x_lambda_525_axial_-350um_to_650um_zstep_10_NA_0.5_fml_393.3um_ftl_265mm_OSR_5_ds_3\psf_zmat';

preprocess_param.num_rawdata = 3000;


%% Parameter setting
estimate_patch_size = 400;
%% Realign parameter
preprocess_param.Nnum = 15; % 15 x 15 pixels behind each microlen. This parameter do not have to change.  
preprocess_param.Nshift = 3; % scanning for Nshift x Nshift times (all choices: 3: 3 x 3; 5: 5 x 5; 13: 13 x 13)
preprocess_param.start_frame = 0; % the number of start frame (the first number is 0, from 0 to N-1)
preprocess_param.frame_interval = 1; % interval of frame (1 by default: process frame by frame, no jump frame)
% Here we choose the center of rawdata and determine ROI, we strongly
% recommend choose the center manually.
preprocess_param.auto_center_mode = 0; % find the center coordinate of x and y automatically when it is 1 and it will disable center_X and center_Y, otherwise 0
preprocess_param.conf_name = './utility/realign/3x3.conf.sk.png'; % configuration file for scanning which is corresponding to Nshift
preprocess_param.center_X = 3999; % the center coordinate x of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.center_Y = 2997; % the center coordinate y of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.Nx = 260; % take half number of microlens in x direction (total number: Nx * 2 + 1)
preprocess_param.Ny = 190; % take half number of microlens in y direction (total number: Nx * 2 + 1)

preprocess_param.group_mode = 1; % Mode of realign between different frame of rawdata. 0: jump mode (group1: 1-9, group2: 10-18,...); 1: slide window(group1: 1-9, group2: 2-10,...)
preprocess_param.slight_resize = 0.9991; % slight resize raw data in realign function (1 by default)
preprocess_param.slight_rotation = 0; % slight rotate raw data in realign function (0 by default) Note that we do not recommend resize and rotate in realign module.

preprocess_param.rotation =  0; % rotate raw data clockwise (all choice: 0, 90, 180, 270)
preprocess_param.upsampling_resize = 0; % 1 means resize WDF to 13 x 13, otherwise no resize; 
preprocess_param.realign_mode = 'LZ'; % realignMode for different scanning sequence (all choices: 'LZ': light path scanning (in RUSH3D). 'ZGX': stage scanning in the opposite direction from 'LZ')
preprocess_param.bright_scale_normalize = 0; % 1 for normalizing bright scale for each scanning period (for lym) 0 no;
preprocess_param.skip_zero_frame = 0; % 1 for skip zero frame, 0 no;
preprocess_param.view_range = 5;
%% Multiscale detrending parameter
debgrecon_param.mul_ds_psf = 30; % down sample rate of psf is 12
debgrecon_param.PSF_broader = 1; % cut psf 12 pixel for each side (four sides)
debgrecon_param.psf_end = 101; % psf axial number
debgrecon_param.zstep = 80; % dz (um)
debgrecon_param.zvrange = [-200,200]; % valid
debgrecon_param.zwrange = [-1000,1000]; % whole
debgrecon_param.view_range = 5;
debgrecon_param.writemode = 0; % 0: no output, 2: output last iteration, 1: ouput all result
debgrecon_param.dispmode = 0; % 0: no disp, 1: disp
debgrecon_param.Nbx = 1;
debgrecon_param.Nby = 1; % Block apart 5 x 5 pieces when DAO
debgrecon_param.maxIter = 1 ; % Max iteration times: 2 or 3 is enough
debgrecon_param.AOstar = 1; % 1 for DAO; 0 for no DAO
debgrecon_param.defocus = 1; % 1 for defocus, 0 for no defocus
debgrecon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
debgrecon_param.margin = 9; % margin overlap
debgrecon_param.mchannel = 1; % number of channel
%% Registration and std parameter
regstd_param.reg_gSig = 2;
regstd_param.reg_bin_width = 200;
regstd_param.rankdetrending = 0;
regstd_param.maxIter = 10;
%% Reconstruction parameter
recon_param.PSF_broader = 50; %% Cut PSF each side for 276 pixels;
recon_param.ds = 3; %% Downsampling rate for PSF
recon_param.psf_end = 101;
recon_param.multi_flag = 0;

recon_param.writemode = 1;
recon_param.dispmode = 1;
recon_param.Nbx = 2;
recon_param.Nby = 2; % Block apart 1 x 1 pieces when DAO
recon_param.num_block = 1 ; % Axial block for 10 when forword propagate
recon_param.maxIter = 3 ; % Max iteration times: 2 or 3 is enough
recon_param.view_range = 5; % About 25 Angle views within the iteration
recon_param.AOstar = 1; % 1 for DAO; 0 for no DAO
recon_param.defocus = 0; % 1 for defocus, 0 for no defocus
recon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
recon_param.margin = 9; % margin overlap

recon_param.patchx = 2;
recon_param.patchy = 2;
recon_param.side = 50;

% this three parameter is not used in reconstruction but use in calcium extraction
recon_param.pixel_size = 4.6e-6/5.761*recon_param.ds; 
recon_param.per_slice_depth = 10e-6; % depth
recon_param.estimate_patch_size = estimate_patch_size;
%% View configuration
view_array = view_config(preprocess_param);
preprocess_param.view_array = view_array;

%% load STD and WDF
% std = loadtiff('D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test1\std_rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_2106221446340.tif');
% std = reshape(std,[size(std,1),size(std,2),preprocess_param.Nnum,preprocess_param.Nnum]);
% std = permute(std,[1,2,4,3]);
% std_vi = zeros(size(std,1),size(std,2),length(view_array));
% for v= 1:length(view_array)
%     std_vi(:,:,v) = std(:,:,view_array{v}(1),view_array{v}(2));
% end
% saveastiff(uint16(std_vi),'D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test2\std_rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_2106221446340.tif');
% std_WDF = loadtiff(sprintf('%s\\std_rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_2106221446340.tif',outdir));
first_wdf = double(loadtiff(sprintf('%s\\realign\\realign_rasgrf-ai148d_5xTMP_whisker3_3x3_50.0ms_Full_Hardware_LaserCount1_210928213137_No0.tif',outdir)));
savegroup = 3;
std_volume = loadtiff(sprintf('%s\\std_recon_seperate\\whole_patch.tiff',outdir));
%%
reg_savepath = 'D:\RUSH3Dproject\RUSH3Dresult\20210928\rasgrf-ai148d_5xTMP_whisker3\test\reg_path';
preprocess_param.valid_frame_num = 2992;
valid_frame_num = preprocess_param.valid_frame_num;
% load psf
disp('---------------------------load psf---------------------------')
recon_param.Nshift = preprocess_param.Nshift; %% Scanning Times = 3 x 3
recon_param.Nnum = preprocess_param.Nnum; %% 15 x 15 pixels behind each MicroLen
recon_param.psf_layer_position = 1:1: recon_param.psf_end;
[psf, recon_param] = psf_load_module(recon_param,recon_psfpath,view_array);
recon_param.upsampling = recon_param.Nnum/recon_param.Nshift/recon_param.ds;
recon_param.patch_side = ceil(recon_param.side *  recon_param.upsampling);


std_volume = double(std_volume);
std_volume = std_volume  / max(std_volume(:));
first_wdf = double(first_wdf);
first_wdf = first_wdf / max(first_wdf(:));
%% start processing patch
NN_ds = preprocess_param.Nnum / preprocess_param.Nshift;
% ------------------------z range-------------------------------
z_energy = squeeze(mean(mean(std_volume,1),2));
start_ind = 13;
end_ind = 98;

std_volume_cut = std_volume(:, :, start_ind : end_ind);
std_volume_cut = std_volume_cut / max(std_volume_cut(:));
volum_size_cut = [size(std_volume_cut, 1), size(std_volume_cut, 2), size(std_volume_cut, 3)];

% --------------------segmentation parameters--------------------
seed_param.outdir = outdir;
seed_param.Nnum = preprocess_param.Nnum;
seed_param.per_slice_depth = recon_param.per_slice_depth;
seed_param.estimate_patch_size = estimate_patch_size;
seed_param.pixel_size = recon_param.pixel_size; % lateral pixel size
seed_param.NN_ds = NN_ds;
seed_param.wdf_ds = recon_param.ds;
seed_param.vol_ds = recon_param.ds;
seed_param.neuron_number = 150;
seed_param.neuron_lateral_size = 6;
seed_param.local_constrast_th = 1.5;
seed_param.optical_psf_ratio = seed_param.per_slice_depth / seed_param.pixel_size;
seed_param.overlap_ratio = 0.5;
seed_param.boundary_thickness = 2; % for plot
seed_param.discard_component_save = false;
seed_param.volume_threshold = [5 * seed_param.neuron_lateral_size.^2, 500* seed_param.neuron_lateral_size.^2];

% --------------------seed generation--------------------
shell_radius = ceil(2* seed_param.neuron_lateral_size);

% --------------------iteration related--------------------
bg_iter=10;
max_demixing_round = 3;
maxIter_NMF = 10;

global_center = [];
global_trace = [];
global_seg = [];

% --------------------background rejection--------------------
oasis_lambda = 0.02;
oasis_g = 0.9;
lambda_l0 = 0.01;
frames_step = 1;
%% vessel segmentation
cv_wdf = squeeze(first_wdf(:,:,ceil(size(first_wdf,3)/2)));
response_stack_segm = vessel_segmentation(cv_wdf, outdir);
response_stack_segm = response_stack_segm(recon_param.side+1:end-recon_param.side,recon_param.side+1:end-recon_param.side);
response_stack_segm_ds = imresize(response_stack_segm,1/seed_param.wdf_ds);
%% patch preparation
% --------------------patch preparation--------------------

wdf_cut_side = first_wdf(recon_param.side+1:end-recon_param.side,recon_param.side+1:end-recon_param.side,:);
wdf_cut_side_ds = imresize(wdf_cut_side,1/seed_param.wdf_ds);
[patch_info_array] = determine_patch_size_LFM(size(std_volume, 1), size(std_volume, 2),...
                    size(wdf_cut_side_ds,1), size(wdf_cut_side_ds, 2), seed_param);
%%
for global_patch_id = 1 :1: 35 % for lateral patches
    
    fprintf('patch %d\n', global_patch_id);
    curr_outdir = sprintf('%s\\patch_%d', outdir, global_patch_id);
    mkdir(curr_outdir);
    %% volume preparation
    curr_patch_info = patch_info_array{global_patch_id};
    patch_volume = std_volume_cut(curr_patch_info.location(1, 1) : curr_patch_info.location(2, 1), ...
        curr_patch_info.location(1, 2) : curr_patch_info.location(2, 2), ...
        :);
    patch_volum_size_cut = [size(patch_volume, 1), size(patch_volume, 2), size(patch_volume, 3)];
    
    patch_vessel_seg = response_stack_segm_ds(curr_patch_info.wdf_loc(1, 1) : curr_patch_info.wdf_loc(2, 1), ...
        curr_patch_info.wdf_loc(1, 2) : curr_patch_info.wdf_loc(2, 2), :);
    patch_wdf_size_cut = [size(patch_vessel_seg,1),size(patch_vessel_seg,2)];
    %% neuron segmentation generation module
    % ~ 2h
	center_array = [];
    disp('--------------------------Neuron segmentation--------------------------')
    % only keep central ones
    curr_seed_param = seed_param;
    curr_seed_param.outdir = curr_outdir;

    valid_seg = segmentation_module(patch_volume, patch_vessel_seg, curr_seed_param);
                                %%
    % determine if it is an empty patch
    if find(~cellfun(@isempty,valid_seg))
        % calcualte component center
        for i = 1 : size(valid_seg, 1)
            center_array(i, :) = mean(valid_seg{i, 2}, 1);
        end

        %% seed module
        % ~ 4h
        % generate seed template for iterations
        disp('--------------------------Iteration preparation--------------------------')
        % define wigners that is usefull
        specify_wigner = zeros(5,1);
        specify_wigner(1) = find(cellfun(@(x)all(x(:)==[ceil(recon_param.Nnum / 2); ceil(recon_param.Nnum / 2)]),view_array));
        for i = 1 : 4
            [u, v] = ind2sub([2, 2], i);
            buf = [ceil(recon_param.Nnum / 2) + (-1)^u * ceil(recon_param.Nnum / 5),...
                   ceil(recon_param.Nnum / 2) + (-1)^v * ceil(recon_param.Nnum / 5)];
            specify_wigner(i+1) = find(cellfun(@(x)all(x(:)==buf(:)),view_array));
        end
        

        % generate 2d initialized seed
        [S_init, S_shell_init,S_mask_init, S_shell_mask_init, valid_seg] = seed_generation_module(psf(:, :, :, start_ind:end_ind), valid_seg, ...
                                                                        patch_volum_size_cut, patch_wdf_size_cut, ...
                                                                        specify_wigner, shell_radius, curr_seed_param);

       

        seg_filtered = valid_seg;
        center_array_filtered  = center_array;

    %% apply patches shifts
        reg_center = center_array_filtered;
        reg_seg = seg_filtered;
        for k = 1 : size(seg_filtered, 1)
            % note each seg has multiple small patches
            reg_seg{k, 2}(:, 1) = seg_filtered{k, 2}(:, 1) + curr_patch_info.location(1, 1);
            reg_seg{k, 2}(:, 2) = seg_filtered{k, 2}(:, 2) + curr_patch_info.location(1, 2);
            
            reg_center(:, 1) = center_array_filtered(:, 1) + curr_patch_info.location(1, 1);
            reg_center(:, 2) = center_array_filtered(:, 2) + curr_patch_info.location(1, 2);
        end
        
        
        
        save(fullfile(curr_outdir, ['raw_registered.mat']),...
            'reg_center', 'reg_seg', '-v7.3');
    else
       continue
    end
end