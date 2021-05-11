clc, clear
close all

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.

%  last update: 5/11/2021. YZ

addpath(genpath('utility'))
outdir = 'out\\crystal_skull_std_dense';
mkdir(outdir)
%%
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

%% data loader
% PSF Path
psfpath = 'Z:/PSF/psf_zww/20200125_genepsf_3.1746x_sim_neg400T400_dz4_15Nnum_OSR5';


% Wigner Path (i.e. raw path)
wdfdata_folder1 = 'D:/RichardW2021/RUSH3DResult/Preprocess/';
wdfdata_folder2 = rawdata_folder2;
wdfdata_folder12 = strcat(wdfdata_folder1, wdfdata_folder2);
wdfdata_folder3 = strcat(wdfdata_folder12, rawdatafilename_perfix,'/');
wdfdata_folder4 = strcat(wdfdata_folder3, 'WDF/');

wdfname = [wdfdata_folder4,'rank1__No0.tif'];
%% parameters
%  reconstruction parameters
param.M = 3.1746; %% Magnification = 3.17
param.Nshift = 3; %% Scanning Times = 3 x 3
param.Nnum = 15; %% 15 x 15 pixels behind each MicroLen
param.PSF_broader = 276; %% Cut PSF each side for 276 pixels;
param.PSF_size = 5; %%£¨x,y,u,v,z£©---->(x,y,u,z) no meaning
param.Block_num = 0; %% no meaning


% ------------------------- Reconstruction parameter ------------------------- % 
param.Angle_range = 25; %% About 25 Angle views within the iteration
param.AOstar = 1; %% 1 for DAO; 0 for no DAO 
param.maxIter = 3; %% Max iteration times: 2 or 3 is enough
param.rotWDF = 0; %% no meaning
param.DownSampling_for_PSF = 1; %% Downsampling rate for PSF
param.UpSamplingRate_for_WDF = Nnum/Nshift/DownSampling_for_PSF; %% UpsamlpingRate for WDF, which is corresponding to the PSF Size.
param.psf_layer = [1:2:51,52:1:151,153:2:201]; % PSF to be loaded []

% ------------------------- DAO parameter ---------------------------
param.defocus = 1; %% 1 for defocus, 0 for no defocus
param.Nbx = 5; param.Nby = 5; %% Block apart 5 x 5 pieces when DAO 
param.num_block= 18 ;    %% Axial block for 10 when forword propagate 
param.threshhold = 25;  %% Shift should not be allowed to exceed [-25,25]
param.sidelobe=9; %% Each 

param.estimate_patch_size = 700;

param.pixel_size = 1.2e-6; % lateral pixel size
param.per_slice_depth = 4e-6; % axial slice depth range

% multi-scale depth calculation
psf_layer_position = [1 : 2 : 51, 52 : 1 : 151, 153 : 2 : 201];

%% realign 
%  require: per-image realignment
realign_param.Nshift = 3;
realign_param.Nnum = 15;

realign_param.large_cycle = 24; %% 24 x 40 Images will be loaded;
realign_param.small_cycle = 40; %% one stack contains 40 LF Images;
realign_param.maxIter = 10; %% Maximum times for Iterate rank normal one

summary_image = preprocessing_module(realign_param, data_path, wdf_path);


%% summary image generation
% current we use rank-1 backgropund detrending + std


%% reconstruction module
curr_volume = reconstruction_module(param,wdfname);
curr_volume = double(curr_volume);
curr_volume = curr_volume  / max(curr_volume(:));

%% neuron segmentation generation module
[global_A_in, global_center] = segmentation_module(curr_volume, seed_param);


%% seed module
% generate seed template for iterations
[seed_array, seed_array_mask] = seed_generation_module(valid_seg, target_wigner);

% note the seed is stored in 
%% background rejection 
%  subtract background informations in the target volume size
%  reconstruct the side view and project to other views?

processed_movie = background_removal(input_movie);


%
%% main iteration module
%  the iteration is based on wigner
max_demixing_round = 5;

for i = 1 : max_demixing_round
	[S, S_bg, bg_spatial] = update_spatial_lasso_with_bg(sensor_movie, S, S_bg, T, ...
            T_bg, bias, bg_spatial, bg_temporal, S_mask_init,S_bg_mask_init, movie_size, maxIter_NMF);        

	[T, T_bg, bg_temporal] = update_temporal_oasis_with_bg(sensor_movie, ...
            S, S_bg, T, T_bg, bias, bg_spatial, bg_temporal, S_mask_init, S_bg_mask_init, movie_size, maxIter_NMF);     
    
end


%% output


%%
figure, subplot(1, 3, 1), spatial_comp_render(obj.A, global_A_in, [movie_size_h, movie_size_w])
subplot(1, 3, 2), obj = CNMFE_show_contours(obj, 0.99); % larger parameter, larger countor
subplot(1, 3, 3), temporal_trace_render(obj.C_raw(1 : 100, :), 'k'), sgtitle('after more runs')
saveas(gca, sprintf('%s\\after_more_runs.png', outdir))


%% final plot: 3D neuron distributions in a video
% with AZ changed
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
	plot_3D_distribution(valid_seg_global_filt_final, [size_h, size_w], [1, 2] * pixel_size, ...
                        param.psf_layer * per_slice_depth, 10) 
	view(az_angle(i), el)
    
	temp = getframe(gcf);
%     temp = imresize(temp.cdata, [400, 600]);
    avi_file.writeVideo(temp);
   
end
avi_file.close();   

%% final plot: temporal activity
tmp_C = zscore(obj.C, 0, 2);
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
%%
figure('position', [1, 1, 300, 600]),
temporal_trace_render_simple(tmp_C(zoom_in_neuron_id(1) : zoom_in_neuron_id(2), ...
                         zoom_in_time(1) : zoom_in_time(2)))
colormap('gray')
caxis([-2, 0])
axis off
pbaspect([1, 2, 1])
%% utility function










