clc;clear;
close all;

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  patched version, for the sake of processing memory and also processing
%  speed.
%  registration module is inserted  after realignment. Dynamic memory calcualtion
%  enabled.



%  last update: 6/29/2021. MW
%  last update: 6/5/2021. YZ
%  last update: 5/23/2021. MW

%% addpath
addpath('background_rejection');
addpath('main demixing');
addpath('../preprocess_module');
addpath('reconstruction_module');
addpath('../registration');
addpath('seed_generation_module');
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
outdir = 'D:\RUSH3Dproject\RUSH3Dresult\20210928\rasgrf-ai148d_5xTMP_whisker3\test1'; % YZ
mkdir(outdir)
outdir0 = 'D:\RUSH3Dproject\RUSH3Dresult\20210928\rasgrf-ai148d_5xTMP_whisker3\test';
first_file_name = 'rasgrf-ai148d_5xTMP_whisker3_3x3_50.0ms_Full_Hardware_LaserCount1_210928213137';
realign_savepath = strcat(outdir0, '\\', 'realign'); 
realign_savename = strcat('realign_', first_file_name); 

realigndata_name_perfix = strcat(realign_savepath, '\\',realign_savename); 

% rawdata dir and file name
% rawdata_path = 'D:\RUSH3Dproject\RUSH3Drawdata\0816\5x_br_tb_0.00_bg_0.00_lfm_uni_N_975_stack'; % raw data path (before rotate, resize and rotate);
rawdata_path = 'D:\RUSH3Dproject\RUSH3Drawdata\20210928\rasgrf-ai148d_5xTMP_whisker3'; % raw data path (before rotate, resize and rotate);



first_file_name = 'rasgrf-ai148d_5xTMP_whisker3_3x3_50.0ms_Full_Hardware_LaserCount1_210928213137'; % in one video or one capture, the name of first stack
input_rawdata_perfix = [rawdata_path, '\\', first_file_name]; % the first file in one capture

% debg psfpath
% debg_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\10x_lambda_488_axial_-600um_to_600um_zstep_20_NA_0.5_fml_875um_ftl_265mm_OSR_5_ds_12\psf_zmat';
debg_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.6509x_axial_-1000um_to_1000um_zstep_20_NA_0.5um_fml_393.3um_ftl_265mm_OSR_5_ds_30\psf_zmat';
recon_psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\5.761x_lambda_525_axial_-300um_to_600um_zstep_10_NA_0.5_fml_393.3um_ftl_265mm_OSR_5_ds_3\psf_zmat';


preprocess_param.video_option = 1; % realign and reconstruction one frame by one frame when 1
preprocess_param.num_rawdata = 100;


%% Parameter setting
estimate_patch_size = 500;
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
debgrecon_param.z_select = 0; % 0 auto select zrange; 1: set manually
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
%% Registration and std parameter
regstd_param.reg_gSig = 2;
regstd_param.reg_bin_width = 200;
regstd_param.rankdetrending = 0;
regstd_param.maxIter = 10;


first_WDF = loadtiff(sprintf('%s_No0.tif',realigndata_name_perfix));

preprocess_param.valid_frame_num = 2992;
valid_frame_num = preprocess_param.valid_frame_num;
debgrecon_param.savegroup = 3;

view_array = view_config(preprocess_param);
preprocess_param.view_array = view_array;

%% Registration and STD
disp('---------------------------Registration---------------------------');

reg_savepath = sprintf('%s\\reg_path', outdir);
if ~exist(reg_savepath,'file')
    mkdir(reg_savepath);
end

std_savepath = sprintf('%s\\std_path', outdir);
if ~exist(std_savepath,'file')
    mkdir(std_savepath);
end

vid_debg_mul_video_savepath = sprintf('%s\\vid_mul_video',outdir0);

t_regstd_start = clock; 

disp('--------------------------Measure shift--------------------------');
tic;
% measure shift
cv_ind = ceil(length(view_array)/2);
centerview_video = [];
shifts = [];
overlap_r = 10;
std_g = zeros(size(first_WDF, 1), size(first_WDF, 2));
t_vone = 0;
for g_id = 1 : debgrecon_param.savegroup
    centerview_video = double(loadtiff(sprintf('%s\\mul_video_view%d_ch1_g%d.tiff',vid_debg_mul_video_savepath, cv_ind, g_id)));
    if g_id > 1
        centerview_video = cat(3,video_overlap,centerview_video);
    end
    [d1,d2,~] = size(centerview_video);
    [~, shifts_seg, bound, option_r] = motion_correction(centerview_video, d1, d2, regstd_param, outdir);
    if g_id > 1
        shifts_seg = shifts_seg(overlap_r+1:end);
        centerview_video = centerview_video(:,:,overlap_r+1:end);
    end
    centerview_video = apply_shifts(centerview_video, shifts_seg, option_r, bound/2, bound/2);
    video_overlap = centerview_video(:,:,end-overlap_r+1:end);
    shifts = cat(1, shifts, shifts_seg);
    t_v = size(centerview_video,3);
    
    % save and std
    saveastiff(im2uint16(centerview_video/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, cv_ind ,g_id));
    
    % calculate std
    maxIter = regstd_param.maxIter;
    if regstd_param.rankdetrending == 1
        [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(centerview_video, [], size(centerview_video,3)), maxIter);
        curr_std_image = compute_std_image(reshape(centerview_video, [], size(centerview_video,3)), ...
            curr_bg_spatial(:), curr_bg_temporal);
    else
        curr_std_image = compute_std_image(reshape(centerview_video, [], size(centerview_video,3)));
    end
    std_g = (t_vone *std_g.^2 + t_v * (reshape(curr_std_image, [size(centerview_video, 1), size(centerview_video, 2)])).^2)./(t_vone+t_v);
    save(sprintf('%s\\std_g%d_v%d.mat',std_savepath, g_id, cv_ind),'std_g','-v7.3');
    t_vone = t_vone + t_v;
    
end
save(sprintf('%s//shifts_wholeprocess.mat',outdir),'shifts','-v7.3');
t_ms = toc;
fprintf('measure shift process done and it takes %.2f secs\n',t_ms);
%%
disp('--------------------------Apply shift--------------------------');


% apply shift for each angle
for v_ind = 2 : length(view_array)
    tic;
    % apply shift
    v = seq_ind(v_ind);
    view_video = [];
    t_vone = 0;
    std_g = zeros(size(first_WDF, 1), size(first_WDF, 2));
    for g_id = 1: debgrecon_param.savegroup
        
        view_video = double(loadtiff(sprintf('%s\\mul_video_view%d_ch1_g%d.tiff', vid_debg_mul_video_savepath, v, g_id)));
        t_v = size(view_video,3);
        shift_seg = shifts(t_vone + 1:t_vone + t_v);
        view_video = apply_shifts(view_video, shifts_seg, option_r, bound/2, bound/2);
        
        saveastiff(im2uint16(view_video/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, v ,g_id));
        
        % calculate std
        maxIter = regstd_param.maxIter;
        if regstd_param.rankdetrending == 1
            [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(view_video, [], size(view_video,3)), maxIter);
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)), ...
                curr_bg_spatial(:), curr_bg_temporal);
        else
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)));
        end
        std_g = (t_vone * std_g.^2 + t_v * (reshape(curr_std_image, [size(view_video, 1), size(view_video, 2)])).^2)./(t_vone+t_v);
        save(sprintf('%s\\std_g%d_v%d.mat',std_savepath, g_id, v),'std_g','-v7.3');
        t_vone = t_vone + t_v;
    end
    t_onereg = toc;
    fprintf('%d in %d view has been registered and std and it takes %.2f secs\n', v, length(view_array), t_onereg);
end
t_regstd = etime(clock,t_regstd_start);
fprintf('resgistration process done and it takes %.2f secs\n',t_regstd); 
%%
% std_concate
std_WDF = zeros(size(first_WDF, 1), size(first_WDF, 2), length(view_array));
for v = 1: length(view_array)
    load(sprintf('%s\\std_g%d_v%d.mat',std_savepath,debgrecon_param.savegroup,v),'std_g');
    std_WDF(:,:,v) = std_g;
end
% save std
saveastiff(im2uint16(std_WDF / max(std_WDF(:))), sprintf('%s\\std_%s%d.tif', outdir, first_file_name,0));