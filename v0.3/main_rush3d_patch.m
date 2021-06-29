clc, clear
close all

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  patched version, for the sake of processing memory and also processing
%  speed.
%  registration module is inserted  after realignment. Dynamic memory calcualtion
%  enabled.

%  last update: 6/29/2021. MW
%  last update: 6/5/2021. YZ
%  last update: 5/23/2021. MW

addpath('background_rejection');
addpath('main demixing');
addpath('preprocess_module');
addpath('reconstruction_module');
addpath('registration');
addpath('seed_generation_module');
addpath(genpath('segmentation_module'));
addpath(genpath('utility'));

outdir = 'D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\\test_patch1';
mkdir(outdir)

rawdata_path = 'D:\RUSH3Dproject\RUSH3Drawdata\0622\rasai148d_1_z287'; % raw data path (before rotate, resize and rotate);
first_file_name = 'rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_210622144634'; % in one video or one capture, the name of first stack
bg_file_name = 'beads_1um_2_bg_3x3_150.0ms_Full_Hardware_LaserCount1_210609000412.0';

input_rawdata_perfix = [rawdata_path, '\\', first_file_name]; % the first file in one capture
bg_rawdata_perfix = [rawdata_path, '\\', bg_file_name]; 

preprocess_param.std_option = 0; % calculating standard deviation of one video and process the std data in the post pipeline when 1
preprocess_param.video_option = 1; % realign and reconstruction one frame by one frame when 1

preprocess_param.subnoise_flag = 0;

preprocess_param.large_cycle = 20;  % total image number: large_cycle * small_cycle
preprocess_param.small_cycle = 40; % In general, one stack has 40 images or less
preprocess_param.Nnum = 15; % 15 x 15 pixels behind each microlen. This parameter do not have to change.  
preprocess_param.Nshift = 3; % scanning for Nshift x Nshift times (all choices: 3: 3 x 3; 5: 5 x 5; 13: 13 x 13)
preprocess_param.mul_ds_psf = 14;
preprocess_param.rotWDF = 0;
%% colormap
% colormap
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

%% Rotate Resize Realign

%  realign will combine the scanned image 
%  the rotation and resize corrected data

disp('---------------------------Realign---------------------------');

realign_savepath = strcat(outdir, '\\', 'realign'); % realign data after realign will be saved in this folder
realign_savename = strcat('realign_', first_file_name); % share the same file name
realigndata_name_perfix = strcat(realign_savepath, '\\',realign_savename); 
if ~exist(realign_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(realign_savepath);
end

preprocess_param.start_frame = 600; % the number of start frame (the first number is 0, from 0 to N-1)
preprocess_param.frame_interval = 1; % interval of frame (1 by default: process frame by frame, no jump frame)
preprocess_param.upsampling_resize = 0;% 1 means resize WDF to 13 x 13, otherwise 0; (has a problem !!!!!)

% Here we choose the center of rawdata and determine ROI, we strongly
% recommend choose the center manually.
preprocess_param.auto_center_mode = 0; % find the center coordinate of x and y automatically when it is 1 and it will disable center_X and center_Y, otherwise 0
preprocess_param.auto_center_frame = preprocess_param.start_frame; % use auto_center_frame to find center automatically under the auto center mode, take the first frame (0 in c++) by default  

preprocess_param.center_X = 4001; % the center coordinate x of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.center_Y = 3000; % the center coordinate y of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.Nx = 255; % take half number of microlens in x direction (total number: Nx * 2 + 1)
preprocess_param.Ny = 196; % take half number of microlens in y direction (total number: Nx * 2 + 1) ( Nx = Ny is strongly recommended) (has a problem !!!!!)

preprocess_param.conf_name = './utility/realign/3x3.conf.sk.png'; % configuration file for scanning which is corresponding to Nshift


preprocess_param.group_mode = 1; % Mode of realign between different frame of rawdata. 0: jump mode (group1: 1-9, group2: 10-18,...); 1: slide window(group1: 1-9, group2: 2-10,...)
if preprocess_param.group_mode == 1
    preprocess_param.group_count = preprocess_param.large_cycle*preprocess_param.small_cycle-preprocess_param.Nshift^2-preprocess_param.start_frame+1; % the number of realigned WDF stacks
else
    preprocess_param.group_count = floor((preprocess_param.large_cycle*preprocess_param.small_cycle-preprocess_param.start_frame)/preprocess_param.Nshift^2);
end


preprocess_param.realign_mode = 'LZ'; % realignMode for different scanning sequence (all choices: 'LZ': light path scanning (in RUSH3D). 'ZGX': stage scanning in the opposite direction from 'LZ')

preprocess_param.rotation =  0; % rotate raw data clockwise (all choice: 0, 90, 180, 270)
preprocess_param.slight_resize = 0.9991; % slight resize raw data in realign function (1 by default)
preprocess_param.slight_rotation = 0; % slight rotate raw data in realign function (0 by default) Note that we do not recommend resize and rotate in realign module.

input_rawdata_name = strcat(input_rawdata_perfix,'.0.tiff');
% start realign
if preprocess_param.video_option == 1
    preprocess_param.centerview = strcat(realigndata_name_perfix,'_cv.tiff');
    first_WDF = realign_module(preprocess_param, input_rawdata_name, realigndata_name_perfix);
end

valid_frame_num = preprocess_param.group_count;

%% load PSF for multiscale
%  reconstruction parameters

psf_param.Nshift = 3; %% Scanning Times = 3 x 3
psf_param.Nnum = 15; %% 15 x 15 pixels behind each MicroLen
psf_param.PSF_broader = 28; %% Cut PSF each side for 276 pixels;
psf_param.multi_flag = 0;

psf_param.PSF_size = 5; %%£¨x,y,u,v,z£©---->(x,y,u,z) no meaning

psf_param.psf_first_Q = 51;
psf_param.psf_second_Q = 151;
psf_param.psf_end = 101;

if psf_param.multi_flag == 1
    psf_param.psf_layer_position = [1 : 2 : psf_param.psf_first_Q, psf_param.psf_first_Q+1 : 1 : psf_param.psf_second_Q, psf_param.psf_second_Q+2 : 2 : psf_param.psf_end];
    psf_param.first_index = find(psf_param.psf_layer_position == psf_param.psf_first_Q);
    psf_param.second_index = find(psf_param.psf_layer_position == psf_param.psf_second_Q);
else
    z_mid_slice_half = 200/20;
    psf_param.psf_layer_position = 1:1: psf_param.psf_end;
    psf_param.first_index = (psf_param.psf_end+1)*0.5-z_mid_slice_half-1;
    psf_param.second_index = (psf_param.psf_end+1)*0.5+z_mid_slice_half;
end

psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-1000um_to_1000um_zstep_20_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5_ds_14\psf_zmat';

[psf, psf_param] = psf_load_module(psf_param,psfpath);

vid_recon_mul_savepath = strcat(outdir, '\\','vid_recon_mul'); 
vid_debg_mul_savepath = strcat(outdir, '\\','vid_mul'); 
if ~exist(vid_recon_mul_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(vid_recon_mul_savepath);
end

if ~exist(vid_debg_mul_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(vid_debg_mul_savepath);
end
%% multiscale detrending
disp('---------------------------multiscale detrending---------------------------');
global bg_ratio
bg_ratio = gpuArray.zeros(psf_param.Nnum,psf_param.Nnum,'single');
for frame_i = 1 : valid_frame_num % time dimension
    tic;

    curr_video = loadtiff(sprintf('%s_No%d.tif', realigndata_name_perfix, frame_i - 1));
    curr_video = single(curr_video);

    % multiscale sub ground noise
    curr_video = sub_mul_bg_module(psf,psf_param,preprocess_param,curr_video,...
        vid_recon_mul_savepath,vid_debg_mul_savepath,first_file_name,frame_i-1);
    
    tt = toc;
    fprintf('%d frame is done, take %.2fs\n',frame_i,tt);
end
clear global bg_ratio;
%% registration and std and multiscale substract ground
disp('---------------------------Registration---------------------------');
reg_path = sprintf('%s\\reg_path', outdir);
mkdir(reg_path)
reg_gSig = 5;
reg_bin_width = 200;
single_img_size = 0.5 / 225 * 80; % in single, gigabite
memory_limt = 0.8e3; % in GB, depent on different PCs
% based on the ram size, decide the cutting
estimate_view_num = memory_limt / ((size(first_WDF, 1) * size(first_WDF, 2)) / 6e3 / 8e3 * valid_frame_num * single_img_size );
estimate_view_num = floor(estimate_view_num);

% registration in a spiral manner
[i_index, j_index] = gen_spiral_center(preprocess_param.Nnum);
i_index = i_index(end : -1 : 1);
j_index = j_index(end : -1 : 1);
reg_group = ceil(preprocess_param.Nnum * preprocess_param.Nnum / estimate_view_num);

% std
maxIter = 10;
max_value = 0;
count_ind = 1;


%%
std_WDF = zeros(size(first_WDF, 1), size(first_WDF, 2), preprocess_param.Nnum^2);
for reg_ind = 1 : reg_group
    curr_view_group_i = i_index((reg_ind - 1) * estimate_view_num + 1: min(reg_ind * estimate_view_num, preprocess_param.Nnum^2));
    curr_view_group_j = j_index((reg_ind - 1) * estimate_view_num + 1: min(reg_ind * estimate_view_num, preprocess_param.Nnum^2));
    
    % read each of the realigned movie
    processed_video = zeros(size(first_WDF, 1), size(first_WDF, 2), valid_frame_num, length(curr_view_group_i), 'single');
    for frame_i = 1 : valid_frame_num % time dimension
        if mod(frame_i, 10) == 0
           fprintf('%d in %d frames collected \n', frame_i, valid_frame_num ) 
        end
        name_perfix = strcat(vid_debg_mul_savepath,'\\wdf_debg_vid');
        curr_video = loadtiff(sprintf('%s_%d.tif',name_perfix, frame_i - 1));
        curr_video = single(curr_video) / 65535;
        
        % multiscale sub ground noise
        curr_video = sub_mul_bg_module(psf,psf_param,preprocess_param,curr_video,...
            vid_recon_mul_savepath,vid_debg_mul_savepath,first_file_name,frame_i-1);
        
        
        
        for view_ind = 1 : length(curr_view_group_i) % number of specified wigner
            buf_ind = sub2ind([preprocess_param.Nnum, preprocess_param.Nnum], ...
                               curr_view_group_i(view_ind), curr_view_group_j(view_ind));
            processed_video(:, :, frame_i, view_ind) = curr_video(:, :, buf_ind);
        end
    end
    
    %% do registration
    % ~ 1h
    if reg_ind == 1
        central_video = processed_video(:, :, :, 1);
        [d1,d2, ~] = size(central_video);
        [~, shifts, bound, option_r] = motion_correction(central_video, d1, d2,reg_gSig, ...
            reg_bin_width, outdir);
        max_value = max(central_video(:)); % note this max_value is not updated
    end
    
    %% apply the registration and calculate std
    for view_ind = 1 : length(curr_view_group_i)
        fprintf('%d in %d view registered \n', count_ind, preprocess_param.Nnum^2)
        curr_video = processed_video(:, :, :, view_ind);
        curr_video = apply_shifts(curr_video, shifts, option_r, bound/2, bound/2);
        
        % save
        buf_ind = sub2ind([preprocess_param.Nnum, preprocess_param.Nnum], ...
            curr_view_group_i(view_ind), curr_view_group_j(view_ind));
        
        saveastiff(imresize(im2uint16(curr_video / max_value / 5),1/4), sprintf('%s\\reg_view_%d.tiff', reg_path, buf_ind));
        
        % calculate std
        [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(curr_video, [], size(curr_video,3)), maxIter);
        curr_std_image = compute_std_image(reshape(curr_video, [], size(curr_video,3)), ...
            curr_bg_spatial(:), curr_bg_temporal);
        std_WDF(:, :, buf_ind) = reshape(curr_std_image, [size(curr_video , 1), size(curr_video , 2)]);
        count_ind = count_ind + 1;
    end
    
end
%
saveastiff(im2uint16(std_WDF / max(std_WDF(:))), sprintf('%s\\std_view.tiff', outdir));
% reshape the std_WDF
std_WDF = reshape(std_WDF, size(std_WDF, 1), size(std_WDF, 2), preprocess_param.Nnum, preprocess_param.Nnum);
std_WDF = permute(std_WDF,[1,2,4,3]);
%% psf parameters
%  reconstruction parameters
disp('---------------------------load psf---------------------------')
psf_param.M = 5.6509; %% Magnification = 3.17
psf_param.Nshift = 3; %% Scanning Times = 3 x 3
psf_param.Nnum = 15; %% 15 x 15 pixels behind each MicroLen
psf_param.PSF_broader = 64; %% Cut PSF each side for 276 pixels;
psf_param.PSF_size = 5; %%£¨x,y,u,v,z£©---->(x,y,u,z) no meaning
psf_param.downsampling_rate = 1; %% Downsampling rate for PSF
psf_param.psf_first_Q = 51;
psf_param.psf_second_Q = 152;
psf_param.psf_end = 101;
psf_param.pixel_size = 4.6e-6/psf_param.M;
psf_param.multi_flag = 0;

if psf_param.multi_flag == 1
    psf_param.psf_layer_position = [1 : 2 : psf_param.psf_first_Q, psf_param.psf_first_Q+1 : 1 : psf_param.psf_second_Q, psf_param.psf_second_Q+2 : 2 : psf_param.psf_end];
else
    psf_param.psf_layer_position = 1:1: psf_param.psf_end;
end

psf_param.per_slice_depth = 4e-6;

psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-200um_to_200um_zstep_4_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5\psf_zmat';

[psf, psf_param] = psf_load_module(psf_param,psfpath);
%%
frame = 0;
wdf_size = size(std_WDF);

preprocess_param.downsampling_rate = 1;
preprocess_param.upsampling = preprocess_param.Nnum / preprocess_param.Nshift / preprocess_param.downsampling_rate;
preprocess_param.rotWDF = 0; %% no meaning

recon_param.writemode = 1;
recon_param.dispmode = 1;
recon_param.Nbx = 1;
recon_param.Nby = 1; % Block apart 1 x 1 pieces when DAO

recon_param.num_block = 1 ; % Axial block for 10 when forword propagate
recon_param.maxIter = 2 ; % Max iteration times: 2 or 3 is enough
recon_param.angle_range = 20; % About 25 Angle views within the iteration


recon_param.AOstar = 1; % 1 for DAO; 0 for no DAO
recon_param.defocus = 1; % 1 for defocus, 0 for no defocus
recon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
recon_param.margin = 9; % margin overlap
recon_param.estimate_patch_size = 700;

recon_param.pixel_size = 1.2e-6; % lateral pixel size
recon_param.per_slice_depth = 4e-6; % axial slice depth range

recon_param.patchx = 5;
recon_param.patchy = 7;
recon_param.side = 29;
recon_param.patch_side = recon_param.side *  preprocess_param.upsampling;

recon_param.patchnum = recon_param.patchx * recon_param.patchy;

index_id_x = round(linspace(recon_param.side+1,wdf_size(1)+1-recon_param.side,recon_param.patchx+1));
index_id_y = round(linspace(recon_param.side+1,wdf_size(2)+1-recon_param.side,recon_param.patchy+1));


x_start = index_id_x(1)-recon_param.side;
x_end   = index_id_x(2)-1+recon_param.side;
y_start = index_id_y(1)-recon_param.side;
y_end   = index_id_y(2)-1+recon_param.side;
std_WDF_sub = std_WDF(x_start:x_end,y_start:y_end,:,:);
std_WDF_up = rot90(imresize(std_WDF_sub, ...
    [floor(size(std_WDF_sub,1)*preprocess_param.upsampling/2)*2+1,...
    floor(size(std_WDF_sub,2)*preprocess_param.upsampling/2)*2+1],'cubic'),2*preprocess_param.rotWDF);

size_xsub_p = size(std_WDF_up,1);
size_ysub_p = size(std_WDF_up,2);
size_xsub = size(std_WDF_up,1)-2*recon_param.patch_side;
size_ysub = size(std_WDF_up,2)-2*recon_param.patch_side;
size_whole_x = size_xsub * recon_param.patchx;
size_whole_y = size_ysub * recon_param.patchy;
std_volume = zeros(size_whole_x,size_whole_y,size(psf,5));

% reconstruction for each patch
for patch_idx = 1: recon_param.patchx
    for patch_idy = 1: recon_param.patchy
        
        patch_id = (patch_idx-1)*recon_param.patchy+patch_idy;
        x_start = index_id_x(patch_idx)-recon_param.side;
        x_end   = index_id_x(patch_idx+1)-1+recon_param.side;
        y_start = index_id_y(patch_idy)-recon_param.side;
        y_end   = index_id_y(patch_idy+1)-1+recon_param.side;
        
        std_WDF_sub = std_WDF(x_start:x_end,y_start:y_end,:,:);
        %% reconstruct with normal deconvolution
        
        disp(['--------reconstruction (',num2str(patch_id),...
                                    '/',num2str(recon_param.patchnum),')---------']);

        
        
        
        std_recon_savepath = strcat(outdir, '\\','std_recon_seperate', '\\','patchid_',num2str(patch_id)); % realign data after realign will be saved in this folder
        std_recon_savename = strcat('std_recon_', first_file_name); % share the same file name
        std_recon_name_perfix = strcat(std_recon_savepath,'\\', std_recon_savename);
        if ~exist(std_recon_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
            mkdir(std_recon_savepath);
        end

        std_WDF_up = rot90(imresize(std_WDF_sub, ...
            [floor(size(std_WDF_sub,1)*preprocess_param.upsampling/2)*2+1,...
            floor(size(std_WDF_sub,2)*preprocess_param.upsampling/2)*2+1],'cubic'),2*preprocess_param.rotWDF);
        std_volume_patch = ones(size(std_WDF_up,1),size(std_WDF_up,2),size(psf,5));
        std_volume_patch = std_volume_patch./sum(std_volume_patch(:)).*sum(std_volume_patch(:))./(size(std_volume_patch,3)*size(std_volume_patch,4));      
        std_volume_patch = reconstruction_module(psf_param, recon_param, std_volume_patch, psf, std_WDF_up, std_recon_name_perfix, frame);
        
        std_volume((patch_idx-1)*size_xsub+1:(patch_idx)*size_xsub,(patch_idy-1)*size_ysub+1:(patch_idy)*size_ysub,:) ...
            = std_volume_patch(recon_param.patch_side+1:size_xsub_p-recon_param.patch_side,...
            recon_param.patch_side+1:size_ysub_p-recon_param.patch_side,:);
        
    end
end
std_recon_savepath = strcat(outdir, '\\','std_recon_seperate');
imwriteTFSK(uint16(std_volume/max(std_volume(:))*65535),strcat(std_recon_savepath,'\\','whole_patch.tiff')); % -- ----- pp
std_volume = double(std_volume);
std_volume = std_volume  / max(std_volume(:));

first_WDF_cut_side = first_WDF(recon_param.side+1:end-recon_param.side,recon_param.side+1:end-recon_param.side,:,:);
std_volume = imresize(std_volume,[size(first_WDF_cut_side,1)*preprocess_param.upsampling,size(first_WDF_cut_side,2)*preprocess_param.upsampling]);
%% start processing patch
% --------------------patch preparation--------------------
estimate_patch_size = round(size(std_volume, 2) / 2);
[patch_info_array]= determine_patch_size_LFM(size(std_volume, 1), size(std_volume, 2),...
                estimate_patch_size, preprocess_param.Nnum);

% --------------------std parameters--------------------
diff_psf_layer_position = diff(psf_param.psf_layer_position);
buf = find(diff_psf_layer_position == 1);
start_ind = buf(1);
end_ind = buf(end);
std_volume_cut = std_volume(:, :, start_ind : end_ind);
std_volume_cut = std_volume_cut / max(std_volume_cut(:));
volum_size_cut = [size(std_volume_cut, 1), size(std_volume_cut, 2), size(std_volume_cut, 3)];


% --------------------segmentation parameters--------------------
seed_param.outdir = outdir;
seed_param.per_slice_depth = psf_param.per_slice_depth;
seed_param.estimate_patch_size = 700;
seed_param.pixel_size = psf_param.pixel_size; % lateral pixel size
seed_param.down_factor = preprocess_param.Nnum / preprocess_param.Nshift;

seed_param.neuron_number = 140;
seed_param.neuron_lateral_size = 13;
seed_param.local_constrast_th = 1.4;
seed_param.optical_psf_ratio = seed_param.per_slice_depth / seed_param.pixel_size;
seed_param.overlap_ratio = 0.5;
seed_param.boundary_thickness = 2; % for plot
seed_param.discard_component_save = false;
seed_param.volume_threshold = [5 * seed_param.neuron_lateral_size.^2, 500* seed_param.neuron_lateral_size.^2];


% --------------------seed generation--------------------
shell_radius = ceil(2* seed_param.neuron_lateral_size);

% --------------------iteration related--------------------
bg_iter=10;
max_demixing_round = 5;
maxIter_NMF = 10;

global_center = [];
global_trace = [];
global_seg = [];

% --------------------background rejection--------------------
oasis_lambda = 0.02;
oasis_g = 0.9;
lambda_l0 = 0.01;
frames_step = 1;

for global_patch_id = 1 : length(patch_info_array) % for lateral patches
    fprintf('patch %d\n', global_patch_id);
    curr_outdir = sprintf('%s\\patch_%d', outdir, global_patch_id);
    mkdir(curr_outdir)
    %% volume preparation
    first_WDF_cut_side = first_WDF(recon_param.side+1:end-recon_param.side,recon_param.side+1:end-recon_param.side,:,:);
    
    curr_patch_info = patch_info_array{global_patch_id};
    patch_volume = std_volume_cut(curr_patch_info.location(1, 1) : curr_patch_info.location(2, 1), ...
                                  curr_patch_info.location(1, 2) : curr_patch_info.location(2, 2), ...
                                  :);
    patch_volum_size_cut = [size(patch_volume, 1), size(patch_volume, 2), size(patch_volume, 3)];
    
    patch_WDF = first_WDF_cut_side(ceil(curr_patch_info.location(1, 1) / ds) : ceil(curr_patch_info.location(2, 1) / ds), ...
                          ceil(curr_patch_info.location(1, 2) / ds) : ceil(curr_patch_info.location(2, 2) / ds), ...
                          :,:);
    %% neuron segmentation generation module
    % ~ 2h
	center_array = [];
    disp('--------------------------Neuron segmentation--------------------------')
    % only keep central ones
    curr_seed_param = seed_param;
    curr_seed_param.outdir = curr_outdir;

    valid_seg = segmentation_module(patch_volume, ...
                                    squeeze(patch_WDF(:, :, ceil(preprocess_param.Nnum / 2), ceil(preprocess_param.Nnum / 2))), ...
                                    curr_seed_param);
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
        specify_wigner = [ceil(psf_param.Nnum / 2), ceil(psf_param.Nnum / 2)];
        for i = 1 : 4
            [u, v] = ind2sub([2, 2], i);
            buf = [ceil(psf_param.Nnum / 2) + (-1)^u * ceil(psf_param.Nnum / 4),...
                   ceil(psf_param.Nnum / 2) + (-1)^v * ceil(psf_param.Nnum / 4)];
            specify_wigner = [specify_wigner ; buf];
        end

        % generate 2d initialized seed
        [S_init, S_shell_init,S_mask_init, S_shell_mask_init] = seed_generation_module(psf(:, :, :, :, start_ind : end_ind), ...
                                                                        psf_param, valid_seg, ...
                                                                        patch_volum_size_cut, ...
                                                                        specify_wigner, shell_radius);

        %% prepare iteration video
        % ~ 2h  
        [processed_video] = iteration_preparation_patch(reg_path, valid_frame_num, specify_wigner, ...
                                                        curr_patch_info, psf_param, S_init);

        % background component initialization
        [bg_spatial_init, bg_temporal_init] = initialize_bg(processed_video, bg_iter);

        % temporal component initialization
        [T_init, T_shell_init] = initialize_T(processed_video, S_mask_init,  S_shell_mask_init);

        %% change S shape
        % ~ 10min
        S = cell_process_A(S_init, size(processed_video)); clear S_init
        S_mask_init= cell_process_A(S_mask_init, size(processed_video)); 
        S_mask_init = S_mask_init > 0;
        % 
        S_bg = cell_process_A(S_shell_init, size(processed_video)) ;  clear S_shell_init
        S_shell_mask_init = cell_process_A(S_shell_mask_init, size(processed_video));
        S_shell_mask_init = S_shell_mask_init > 0;

        save(sprintf('%s\\initialization.mat', curr_outdir), 'S', 'S_bg', 'S_mask_init', ...
                                                        'S_shell_mask_init', 'T_init', 'T_shell_init', '-v7.3')
        %% main iteration module
        disp('--------------------------main NMF--------------------------')
        T = T_init;
        T_bg = T_shell_init;
        bg_spatial =  bg_spatial_init;
        bg_temporal = bg_temporal_init;

        for i = 1 : max_demixing_round
            fprintf('main demixing %d in %d \n', i, max_demixing_round)
            [S, S_bg, bg_spatial] = update_spatial_lasso(processed_video, S, S_bg, T, ...
                    T_bg, bg_spatial, bg_temporal, S_mask_init,S_shell_mask_init, maxIter_NMF);        
        %
            [T, T_bg, bg_temporal] = update_temporal_lasso(processed_video, ...
                    S, S_bg, T, T_bg, bg_spatial, bg_temporal, S_mask_init, S_shell_mask_init, maxIter_NMF);     

        end

        save(sprintf('%s\\final_component.mat', curr_outdir), ...
                'S', 'S_bg', 'T', 'T_bg', 'bg_spatial', 'bg_temporal', '-v7.3')  
        %% further core-shell demxing
        disp('--------------------------background subtraction--------------------------')
        tic

        neuron_trace_mat = zeros(size(T, 1), size(T, 2));
        deconvol_neuron_trace_mat = zeros(size(T, 1), size(T, 2));
        spike_mat = zeros(size(T, 1), size(T, 2));
        coefs_array = zeros(1, size(T, 1));

        if ~isempty(S_bg)
            parfor i = 1 : size(S, 2)      
                center_trace = T(i, :);
                bg_trace = T_bg(i, :);
                [coefs, sub_trace, ca_out, sp_out] = neuropil_coefficient_estimation_greedy(center_trace, ...
                    bg_trace, oasis_lambda, oasis_g, lambda_l0);       
                neuron_trace_mat(i, :) = sub_trace; % there is no change of 
                deconvol_neuron_trace_mat(i, :)  = ca_out;
                spike_mat(i, :) = sp_out;  
                coefs_array(i) = coefs;
            end   
        else
            parfor i = 1 : size(S, 2)   
                center_trace = T(i, :);
                [ca_out, sp_out, ~, ~, ~] = foopsi_oasisAR1(center_trace, oasis_g, oasis_lambda, false,...
                 true, 1, 100);
                neuron_trace_mat(i, :) = center_trace; 
                deconvol_neuron_trace_mat(i, :) = ca_out;
                spike_mat(i, :) = sp_out;     
            end
            coefs_array = [];
        end
        toc
        save(fullfile(curr_outdir, ['after_background_rejection.mat']),...
            'neuron_trace_mat', 'deconvol_neuron_trace_mat', 'spike_mat', 'coefs_array', '-v7.3');
        figure, temporal_trace_render(neuron_trace_mat, 'k')

        %% neuron trace filtering
        disp('--------------------------temporal filtering--------------------------')
        method = 'svm';
        trace_keep_mode = 'sensitive';
        positive_ind = temporal_activity_filtering(neuron_trace_mat, method, trace_keep_mode);
        % update
        T_filtered = neuron_trace_mat(positive_ind, :);
        S_filtered = spike_mat(positive_ind, :);
        T_deconv_filtered = deconvol_neuron_trace_mat(positive_ind, :);
        center_array_filtered  = center_array(positive_ind, :);
        seg_filtered = valid_seg(positive_ind, :);

        save(fullfile(curr_outdir, ['after_classifier_filter.mat']),...
            'T_filtered', 'S_filtered', 'T_deconv_filtered', ...
            'center_array_filtered', 'seg_filtered', '-v7.3');
        
    %% apply patches shifts
        reg_center = center_array_filtered;
        reg_seg = seg_filtered;
        for k = 1 : size(seg_filtered, 1)
            % note each seg has multiple small patches
            reg_seg{k, 2}(:, 1) = valid_seg{k, 2}(:, 1) + curr_patch_info.location(1, 1);
            reg_seg{k, 2}(:, 2) = valid_seg{k, 2}(:, 2) + curr_patch_info.location(1, 2);
            
            reg_center(:, 1) = center_array_filtered(:, 1) + curr_patch_info.location(1, 1);
            reg_center(:, 2) = center_array_filtered(:, 2) + curr_patch_info.location(1, 2);
        end
        
        global_center = [global_center; reg_center];
        global_seg = [global_seg; reg_seg];
        global_trace = [global_trace; T_filtered];
        
        
        save(fullfile(curr_outdir, ['registered.mat']),...
            'reg_center', 'reg_seg', ' T_filtered', '-v7.3');
    else
        continue
    end
    

end
% save global results
save(fullfile(outdir, ['global_output.mat']),...
        'global_center', 'global_seg', 'global_trace', '-v7.3');
%% 3D distributions
figure('Name', 'video', 'position', [100, 100, 400, 300])
plot_3D_distribution(global_seg, [volum_size_cut(1), volum_size_cut(2)], [1, 2] * seed_param.pixel_size, ...
                        seed_param.per_slice_depth * psf_param.psf_layer_position, 10) 
savefig(sprintf('%s\\spatial_distribution_3D.fig', outdir))

%% final plot: 3D neuron distributions in a video
% with AZ changed
%{
avi_filename =[outdir, '\\spatial_video.avi'];
avi_file = VideoWriter(avi_filename);
avi_file.FrameRate = 30;
avi_file.Quality = 100;
% avi_file.CompressionRatio = 1;
% avi_file.VideoCompressionMethod = 'H.264';
avi_file.open();


az_angle = 1 : 720;
el = 19;
figure('Name', 'video', 'position', [100, 100, 400, 300])
for i = 1 : 1 : length(az_angle)
    clf
	plot_3D_distribution(valid_seg, [volum_size(1), volum_size(2)], [1, 2] * pixel_size, ...
                        recon_param.psf_layer * per_slice_depth, 10) 
	view(az_angle(i), el)
    
	temp = getframe(gcf);
%     temp = imresize(temp.cdata, [400, 600]);
    avi_file.writeVideo(temp);
   
end
avi_file.close();   
%}
%% final plot: temporal activity
tmp_C = zscore(global_trace, 0, 2);
figure, imagesc(tmp_C)
pbaspect([1, 2, 1])

% zoom in,
zoom_in_neuron_id = [1000, 1100];
zoom_in_time = [500, 800];
hold on,
rectangle('Position', [zoom_in_time(1), zoom_in_neuron_id(1), ...
                zoom_in_time(2) - zoom_in_time(1), ...
                zoom_in_neuron_id(2) - zoom_in_neuron_id(1)], ...
                'edgecolor', color_scheme_npg(1, :), ...
                'linestyle', '--', ...
                'linewidth', 1)
savefig(sprintf('%s\\temporal_trace_heat.fig', outdir))
% zoom in
figure('position', [1, 1, 300, 600]),
temporal_trace_render_simple(global_trace(zoom_in_neuron_id(1) : zoom_in_neuron_id(2), ...
                         zoom_in_time(1) : zoom_in_time(2)))
colormap('gray')
caxis([-2, 0])
axis off
pbaspect([1, 2, 1])
savefig(sprintf('%s\\temporal_trace.fig', outdir))


%% utility function
function out_mat_A = cell_process_A(A, movie_size) % would be slow since reach the limit
    out_mat_A = zeros(movie_size(1) * movie_size(2) * movie_size(4), size(A, 2));
    for i = 1 : size(A, 1) % number of component
%         i
       curr_A =  cell2mat(A(i, :));
       out_mat_A(:, i) = curr_A(:);
    end
    out_mat_A  = sparse(out_mat_A );
end
function A_cell = extract_A(A_mat, movie_size)
    A_cell = cell(size(A_mat, 1), movie_size(4));
    
	for i = 1 : size(A_mat, 1) % number of component
       for j = 1 : movie_size(4) % number of view
          curr_A = A_mat((j - 1) * movie_size(1) * movie_size(2) + 1 : j * movie_size(1) * movie_size(2), i);
          A_cell{i, j} = sparse(reshape(curr_A, movie_size(1), movie_size(2)));
       end
    end
end