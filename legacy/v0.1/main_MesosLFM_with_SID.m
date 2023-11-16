clc, clear
close all

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  last update: 5/26/2021. YZ
%  last update: 5/23/2021. MW

addpath('background_rejection');
addpath('main demixing');
addpath('preprocess_module');
addpath('reconstruction_module');
addpath('seed_generation_module');
addpath(genpath('segmentation_module'));
addpath(genpath('utility'));

outdir = 'D:\\data\\RUSH3D_GYD\\0309_out';
mkdir(outdir)

rawdata_path = 'D:\data\RUSH3D_GYD\0309'; % raw data path (before rotate, resize and rotate);
first_file_name = 'X1_NORMAL_3x3_25.0ms_HBar_Hardware_LaserCount1_210309202631'; % in one video or one capture, the name of first stack
input_rawdata_perfix = [rawdata_path, '\\', first_file_name]; % the first file in one capture
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


%% psf parameters
%  reconstruction parameters
disp('---------------------------load psf---------------------------')
psf_param.M = 3.1746; %% Magnification = 3.17
psf_param.Nshift = 3; %% Scanning Times = 3 x 3
psf_param.Nnum = 15; %% 15 x 15 pixels behind each MicroLen
psf_param.PSF_broader = 276; %% Cut PSF each side for 276 pixels;
psf_param.PSF_size = 5; %%£¨x,y,u,v,z£©---->(x,y,u,z) no meaning
psf_param.downsampling_rate = 1; %% Downsampling rate for PSF
psf_param.psf_layer_position = [1 : 2 : 51, 52 : 1 : 151, 153 : 2 : 201];

psfpath = 'D:\data\RUSH3D_GYD\20200125_genepsf_3.1746x_sim_neg400T400_dz4_15Nnum_OSR5';

[psf, psf_param] = psf_load_module(psf_param,psfpath);

%% rotate && resize
disp('---------------------------rotate_resize---------------------------')
rr_savepath = strcat(outdir,'\\rotate_resize'); % raw data after rotate and resize will be save as the same format in this folder
rr_savename = strcat('rr_', first_file_name); % share the same file name
rr_rawdata_name = strcat(rr_savepath, '\\', rr_savename); 
if ~exist(rr_savepath,'file')
    mkdir(rr_savepath);
end

preprocess_param.downsampling_rate = 1; % Downsampling rate for PSF
preprocess_param.large_cycle = 24;  % total image number: large_cycle * small_cycle
preprocess_param.small_cycle = 40; % In general, one stack has 40 images or less
preprocess_param.pre_rotate = -0.029; % slight rotate (degree),
preprocess_param.pre_resize = 1; % slight resize (proportion)
preprocess_param.Nnum = 15; % 15 x 15 pixels behind each microlen. This parameter do not have to change.  
preprocess_param.Nshift = 3; % scanning for Nshift x Nshift times (all choices: 3: 3 x 3; 5: 5 x 5; 13: 13 x 13)
preprocess_param.upsampling = preprocess_param.Nnum / preprocess_param.Nshift / preprocess_param.downsampling_rate;
preprocess_param.rotWDF = 0; %% no meaning

[raw_stack_after_rr, total_frame_num] = rotate_resize_module(preprocess_param, input_rawdata_perfix, rr_rawdata_name); % correct raw size

%% std with detrending(optional)

% Note that at least one of them is 1.
preprocess_param.std_option = 1; % calculating standard deviation of one video and process the std data in the post pipeline when 1

preprocess_param.maxIter = 10;

std_savepath = strcat(outdir, '\\', 'std'); % realign data after realign will be saved in this folder
std_savename = strcat('std_', first_file_name); % share the same file name
std_data_name = strcat(std_savepath, '\\',std_savename); 
if ~exist(std_savepath,'file')
    mkdir(std_savepath);
end

if preprocess_param.std_option == 1
    raw_stack_after_std = std_module(preprocess_param, raw_stack_after_rr,std_data_name); % 
end

%% realign raw data and std data 
%  realign will combine the scanned image 
%  the rotation and resize corrected data
disp('---------------------------Realign---------------------------');
rawdata_name = strcat(rr_rawdata_name,'.',num2str(0),'.tiff'); % in one video or one capture, the name of first stack (***.0.tiff) including folder path

realign_savepath = strcat(outdir, '\\', 'realign'); % realign data after realign will be saved in this folder
realign_savename = strcat('realign_', first_file_name); % share the same file name
realigndata_name_perfix = strcat(realign_savepath, '\\',realign_savename); 
if ~exist(realign_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(realign_savepath);
end


std_realign_savepath = strcat(outdir, '\\', 'std_realign');
std_realign_savename = strcat('std_realign_', first_file_name); % share the same file name
std_realigndata_name_perfix = strcat(std_realign_savepath, '\\', std_realign_savename);
if ~exist(std_realign_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(std_realign_savepath);
end

preprocess_param.start_frame = 0; % the number of start frame (the first number is 0, from 0 to N-1)
preprocess_param.frame_interval = 1; % interval of frame (1 by default: process frame by frame, no jump frame)
preprocess_param.upsampling_resize = 0;% 1 means resize WDF to 13 x 13, otherwise 0; (has a problem !!!!!)

% Here we choose the center of rawdata and determine ROI, we strongly
% recommend choose the center manually.
preprocess_param.auto_center_mode = 0; % find the center coordinate of x and y automatically when it is 1 and it will disable center_X and center_Y, otherwise 0
preprocess_param.auto_center_frame = preprocess_param.start_frame; % use auto_center_frame to find center automatically under the auto center mode, take the first frame (0 in c++) by default  
preprocess_param.center_X = 2001; % the center coordinate x of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.center_Y = 2002; % the center coordinate y of Light field data (Note that the coordinate is from 0 to N-1 in c++)
preprocess_param.Nx = 129; % take half number of microlens in x direction (total number: Nx * 2 + 1)
preprocess_param.Ny = 129; % take half number of microlens in y direction (total number: Nx * 2 + 1) ( Nx = Ny is strongly recommended) (has a problem !!!!!)


preprocess_param.conf_name = './utility/realign/3x3.conf.sk.png'; % configuration file for scanning which is corresponding to Nshift


preprocess_param.group_mode = 1; % Mode of realign between different frame of rawdata. 0: jump mode (group1: 1-9, group2: 10-18,...); 1: slide window(group1: 1-9, group2: 2-10,...)
if preprocess_param.group_mode == 1
    preprocess_param.group_count = preprocess_param.large_cycle*preprocess_param.small_cycle-preprocess_param.Nshift^2; % the number of realigned WDF stacks
else
    preprocess_param.group_count = floor(preprocess_param.large_cyclepreprocess_param.*small_cycle/preprocess_param.Nshift^2);
end


preprocess_param.realign_mode = 'LZ'; % realignMode for different scanning sequence (all choices: 'LZ': light path scanning (in RUSH3D). 'ZGX': stage scanning in the opposite direction from 'LZ')
preprocess_param.centerview = 'center.tiff'; % file name of centerview in the same output directory. (If centerview is not needed, use 'None'. If we need centerview only, use *.only.tiff)

preprocess_param.video_option = 1; % realign and reconstruction one frame by one frame when 1

preprocess_param.rotation =  0; % rotate raw data clockwise (all choice: 0, 90, 180, 270)
preprocess_param.slight_resize = 1; % slight resize raw data in realign function (1 by default)
preprocess_param.slight_rotation = 0; % slight rotate raw data in realign function (0 by default) Note that we do not recommend resize and rotate in realign module.

% start realign
if preprocess_param.video_option == 1
    first_WDF = realign_module(preprocess_param, rawdata_name, realigndata_name_perfix);
end

% processing std
if preprocess_param.std_option == 1
    preprocess_param.group_count = 1;
    std_WDF = realign_module(preprocess_param, strcat(std_data_name, '.tiff'), std_realigndata_name_perfix);
end

% note the output of realign will be the number of frameN - shift^2
valid_frame_num = total_frame_num - preprocess_param.Nshift * preprocess_param.Nshift;
movie_size = [size(std_WDF, 1), size(std_WDF, 2), valid_frame_num];

%%
% ~ 2h
disp('---------------------------Reconstruction---------------------------');
gpuDevice(1);
realign_datafilename = strcat(realigndata_name_perfix,'_No');

std_recon_savepath = strcat(outdir, '\\','std_recon'); % realign data after realign will be saved in this folder
std_recon_savename = strcat('std_recon_', first_file_name); % share the same file name
std_recon_name_perfix = strcat(std_recon_savepath, '\\', std_recon_savename); 
if ~exist(std_recon_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(std_recon_savepath);
end

vid_recon_savepath = strcat(outdir, '\\','vid_recon'); % realign data after realign will be saved in this folder
vid_recon_savename = strcat('vid_recon_', first_file_name); % share the same file name
vid_recon_name_perfix = strcat(vid_recon_savepath, '\\', vid_recon_savename); 
if ~exist(vid_recon_savepath,'file') % output file name including folder path (Note that the folder must exist and it will not report an error when missing folder)
    mkdir(vid_recon_savepath);
end

recon_param.angle_range = 25; % About 25 Angle views within the iteration
recon_param.AOstar = 1; % 1 for DAO; 0 for no DAO 
recon_param.maxIter = 3; % Max iteration times: 2 or 3 is enough
recon_param.defocus = 1; % 1 for defocus, 0 for no defocus
recon_param.Nbx = 5; 
recon_param.Nby = 5; % Block apart 5 x 5 pieces when DAO 
recon_param.num_block= 18 ; % Axial block for 10 when forword propagate 
recon_param.threshhold = 25; % Shift should not be allowed to exceed [-25,25]
recon_param.margin = 9; % margin overlap

% reconstruct STD
if preprocess_param.std_option == 1
    frame = 1;
    std_WDF = rot90(imresize(std_WDF, ...
        [floor(size(std_WDF,1)*preprocess_param.upsampling),floor(size(std_WDF,2)*preprocess_param.upsampling)],'cubic'),2*preprocess_param.rotWDF);
    std_volume = ones(size(std_WDF,1),size(std_WDF,2),size(psf,5));
    std_volume = std_volume./sum(std_volume(:)).*sum(std_volume(:))./(size(std_volume,3)*size(std_volume,4));
    std_volume = reconstruction_module(psf_param, recon_param, std_volume, psf, std_WDF, ...
                                        std_recon_savepath, std_recon_name_perfix,frame);
    std_volume = double(std_volume);
    std_volume = std_volume  / max(std_volume(:));
end

% note std_volume is the final output.
volum_size = [size(std_volume, 1), size(std_volume, 2), size(std_volume, 3)];
%% neuron segmentation generation module
% ~ 2h
disp('--------------------------Neuron segmentation--------------------------')
% only keep central ones
diff_psf_layer_position = diff(psf_param.psf_layer_position);
buf = find(diff_psf_layer_position == 1);
start_ind = buf(1);
end_ind = buf(end);
std_volume_cut = std_volume(:, :, start_ind : end_ind);
std_volume_cut = std_volume_cut / max(std_volume_cut(:));
volum_size_cut = [size(std_volume_cut, 1), size(std_volume_cut, 2), size(std_volume_cut, 3)];
seed_param.outdir = outdir;
seed_param.per_slice_depth = recon_param.per_slice_depth;
seed_param.estimate_patch_size = 700;
seed_param.pixel_size = 1.2e-6; % lateral pixel size
seed_param.per_slice_depth = 4e-6; % axial slice depth range ???
seed_param.down_factor = preprocess_param.Nnum / preprocess_param.Nshift;

% segmentation related
seed_param.neuron_number = 140;
seed_param.neuron_lateral_size = 13;
seed_param.local_constrast_th = 1.4;
seed_param.optical_psf_ratio = seed_param.per_slice_depth / seed_param.pixel_size;
seed_param.overlap_ratio = 0.5;
seed_param.boundary_thickness = 2; % for plot
seed_param.discard_component_save = false;
seed_param.volume_threshold = [5 * seed_param.neuron_lateral_size.^2, 500* seed_param.neuron_lateral_size.^2];

valid_seg = segmentation_module(std_volume_cut, ...
                                            first_WDF(:, :, ceil(preprocess_param.Nnum^2 / 2)), ...
                                            seed_param);


%% seed module
% ~ 4h
% generate seed template for iterations
shell_radius = ceil(2* seed_param.neuron_lateral_size);
disp('--------------------------Iteration preparation--------------------------')

% define wigners that is usefull
specify_wigner = [ceil(psf_param.Nnum / 2), ceil(psf_param.Nnum / 2)];
for i = 1 : 4
    [u, v] = ind2sub([2, 2], i);
    buf = [ceil(psf_param.Nnum / 2) + (-1)^u * ceil(psf_param.Nnum / 4),...
           ceil(psf_param.Nnum / 2) + (-1)^v * ceil(psf_param.Nnum / 4)];
    specify_wigner = [specify_wigner ; buf];
end


[S_init, S_shell_init,S_mask_init, S_shell_mask_init] = seed_generation_module(psf(:, :, :, :, start_ind : end_ind), psf_param, valid_seg, ...
                                                                volum_size_cut, specify_wigner, shell_radius);


%% prepare iteration video
% ~ 2h
bg_iter=10;
[processed_video] = iteration_preparation(realigndata_name_perfix, valid_frame_num, specify_wigner);

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

save(sprintf('%s\\initialization.mat', outdir), 'S', 'S_bg', 'S_mask_init', ...
                                                'S_shell_mask_init', 'T_init', 'T_shell_init', '-v7.3')
%% main iteration module
disp('--------------------------main NMF--------------------------')
%  the iteration is based on wigner
max_demixing_round = 5;
maxIter_NMF = 10;

% main demixing
T = T_init;
T_bg = T_shell_init;
bg_spatial =  bg_spatial_init; % 
bg_temporal = bg_temporal_init;


% S = S(:, 1 : 20);
% S_mask_init = S_mask_init(:, 1 : 20);
% S_bg = S_bg(:, 1 : 20);
% S_shell_mask_init = S_shell_mask_init(:, 1 : 20);
% T = T(1 : 20, :);
% T_bg = T_bg(1 : 20, :);

for i = 1 : max_demixing_round
    fprintf('main demixing %d in %d \n', i, max_demixing_round)
	[S, S_bg, bg_spatial] = update_spatial_lasso(processed_video, S, S_bg, T, ...
            T_bg, bg_spatial, bg_temporal, S_mask_init,S_shell_mask_init, maxIter_NMF);        
%
	[T, T_bg, bg_temporal] = update_temporal_lasso(processed_video, ...
            S, S_bg, T, T_bg, bg_spatial, bg_temporal, S_mask_init, S_shell_mask_init, maxIter_NMF);     
    
end

save(sprintf('%s\\final_component.mat', outdir), ...
        'S', 'S_bg', 'T', 'T_bg', 'bg_spatial', 'bg_temporal', '-v7.3')  
%% further core-shell demxing
disp('--------------------------background subtraction--------------------------')
tic
oasis_lambda = 0.02;
oasis_g = 0.9;
lambda_l0 = 0.01;
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
save(fullfile(outdir, ['after_background_rejection.mat']),...
    'neuron_trace_mat', 'deconvol_neuron_trace_mat', 'spike_mat', 'coefs_array', '-v7.3');
figure, temporal_trace_render(neuron_trace_mat, 'k')

%% 3D distributions
figure('Name', 'video', 'position', [100, 100, 400, 300])
plot_3D_distribution(valid_seg, [volum_size(1), volum_size(2)], [1, 2] * seed_param.pixel_size, ...
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
tmp_C = zscore(T, 0, 2);
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
% zoom in
figure('position', [1, 1, 300, 600]),
temporal_trace_render_simple(T(zoom_in_neuron_id(1) : zoom_in_neuron_id(2), ...
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