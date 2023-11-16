clc;clear;
close all;

%% this file is the main file for the meso-sLFM calcium data processing.
%  this pipeline contains from reading to calcium sensing.
%  patched version, for the sake of processing memory and also processing
%  speed.
%  registration module is inserted  after realignment. Dynamic memory calcualtion
%  enabled.


%  last update: 10/19/2021. MW
%  last update: 6/29/2021. MW
%  last update: 6/5/2021. YZ
%  last update: 5/23/2021. MW

%%  Path and PSF path
% addpath
addpath('../background_rejection');
addpath('../main demixing');
addpath('../preprocess_module');
addpath('../reconstruction_module');
addpath('../registration');
addpath('../seed_generation_module');
addpath(genpath('../segmentation_module'));
addpath(genpath('../utility'));

%% load mat

% load param
load('..\\param.mat','main_param'); %%%
outdir = main_param.outdir;

% patch info
load(sprintf('%s\\patch_info.mat',outdir),'patch_info_array');
patch_id = 12; %1-12
patch_info = patch_info_array{patch_id};
h = patch_info.size_ov(1);
w = patch_info.size_ov(2);
tl = patch_info.location_ov(1,:);
rd = patch_info.location_ov(2,:);

h_ds = patch_info.size_ds(1);
w_ds = patch_info.size_ds(2);
tl_ds = patch_info.location_ds(1,:);
rd_ds = patch_info.location_ds(2,:);

outdir_p = sprintf('%s\\patch_%d',outdir,patch_id);

%% multi scale debg for STD Nshift = 3
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'realign_param','psf_param','debg_param');
% multi scale debg param
realigndata_name_perfix = realign_param.realigndata_name_perfix;
debg_param.Nshift = realign_param.Nshift; 
debg_param.Nnum = realign_param.Nnum; 
debg_param.psf_layer_position = 1 : 4 : debg_param.psf_end;
view_array = view_config(main_param);
debg_param.view_array = view_array;
first_file_name = main_param.first_file_name;
maxframe = debg_param.maxframe;

% group in time lapse
savegroup = ceil(realign_param.valid_frame_num/maxframe);


% load psf
debg_psfpath = psf_param.debg_psfpath;
[depsf, debg_param] = psf_load_module(debg_param,debg_psfpath,view_array);
fprintf('load low resolution psf process done\n');

% multiscale detrending start
disp('---------------------------multiscale detrending---------------------------');

% output file name
recon_mul_savepath = sprintf('%s\\recon_mul',outdir_p); 
debg_mul_video_savepath = sprintf('%s\\mul_video',outdir_p);
if ~exist(recon_mul_savepath,'file') 
    mkdir(recon_mul_savepath);
end
if ~exist(debg_mul_video_savepath,'file') 
    mkdir(debg_mul_video_savepath);
end

% bg ratio and z_range
view_video = zeros(h, w, length(view_array));
for v = 1: length(view_array)
    view_video(:,:,v) = double(imread(sprintf('%s_No%d.tif', realigndata_name_perfix, 0),...
        sub2ind([main_param.Nnum,main_param.Nnum],view_array{v}(2),view_array{v}(1)),...
        'PixelRegion',{[tl(1),rd(1)],[tl(2),rd(2)]}));
end
debg_param = bg_param_def(depsf, debg_param,view_video, recon_mul_savepath, first_file_name);


% start de background
t_debg_start = clock;
[~, ~, seq_ind] = Spiral_circle(main_param.Nnum,main_param.view_range);
seq_ind = seq_ind(end:-1:1);

count = 0;
bg_ratio = gpuArray(single(debg_param.bg_ratio));
% debg frame by frame
for frame_i = 1 : realign_param.valid_frame_num
    view_video = zeros(h, w, length(view_array));
    % load image
    for v = 1: length(view_array)
        view_video(:,:,v) = double(imread(sprintf('%s_No%d.tif', realigndata_name_perfix, frame_i - 1),...
            sub2ind([main_param.Nnum,main_param.Nnum],view_array{v}(2),view_array{v}(1)),...
            'PixelRegion',{[tl(1),rd(1)],[tl(2),rd(2)]}));
    end
    % multiscale sub ground noise
    view_video = sub_mul_bg_module(depsf,debg_param,view_video,...
        recon_mul_savepath,first_file_name,frame_i-1,bg_ratio);
    
    % save
    if mod(frame_i-1,maxframe) == 0
        count = count + 1;
        for v_ind = 1:length(view_array)
            v = seq_ind(v_ind);
            imwrite(im2uint16(view_video(:,:,v)/65535),strcat(debg_mul_video_savepath,'\\',...
                'mul_video_view',num2str(v),'_g',num2str(count),'.tiff'));
        end
    else
        for v_ind = 1:length(view_array)
            v = seq_ind(v_ind);
            imwrite(im2uint16(view_video(:,:,v)/65535),strcat(debg_mul_video_savepath,'\\',...
                'mul_video_view',num2str(v),'_g',num2str(count),'.tiff'),'WriteMode','append');
        end
    end
    tt = toc;
    fprintf('%d frame is done, take %.2fs\n',frame_i,tt);
end

t_debg = etime(clock,t_debg_start);
fprintf('Debg process done and it takes %.2f secs\n',t_debg);
clear bg_ratio;
clear depsf;

% save parameters
debg_param.savegroup = savegroup;
debg_param.debg_mul_video_savepath = debg_mul_video_savepath;
save(sprintf('%s\\param_debg.mat',outdir_p), 'debg_param', '-v7.3');
%% Registration and STD
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'reg_param');
load(sprintf('%s\\param_debg.mat',outdir_p),'debg_param');
debg_mul_video_savepath = debg_param.debg_mul_video_savepath;
view_array = view_config(main_param);
disp('---------------------------Registration for Nshift = 3---------------------------');
% save path
reg_savepath = sprintf('%s\\reg_path', outdir_p);
std_savepath = sprintf('%s\\std_path', outdir_p);

if ~exist(reg_savepath,'file')
    mkdir(reg_savepath);
end
if ~exist(std_savepath,'file')
    mkdir(std_savepath);
end


disp('--------------------------Measure shift for Nshift = 3--------------------------');
tic;
% measure shift
cv_ind = ceil(length(view_array)/2);
centerview_video = [];

t_vone = 0;
std_g = zeros(h, w);

for g_id = 1 : debg_param.savegroup
    load_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff',debg_mul_video_savepath, cv_ind, g_id)));
    centerview_video = cat(3,centerview_video,load_video);
end

[d1,d2,~] = size(centerview_video);
[~, shifts, bound, option_r] = motion_correction(centerview_video, d1, d2, reg_param, outdir_p);

% apply shift to center view
maxIter = reg_param.maxIter;
for g_id = 1:debg_param.savegroup
    centerview_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff',debg_mul_video_savepath, cv_ind, g_id)));
    t_v = size(centerview_video,3);
    shifts_seg = shifts(t_vone + 1 : t_vone + t_v);
    centerview_video = apply_shifts(centerview_video, shifts_seg, option_r, bound/2, bound/2);
    % save
    saveastiff(im2uint16(centerview_video/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, cv_ind ,g_id));
    
    % calculate std
    if reg_param.rankdetrending == 1
        [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(centerview_video, [], size(centerview_video,3)), maxIter);
        curr_std_image = compute_std_image(reshape(centerview_video, [], size(centerview_video,3)), ...
            curr_bg_spatial(:), curr_bg_temporal);
    else
        curr_std_image = compute_std_image(reshape(centerview_video, [], size(centerview_video,3)));
    end
    std_g = sqrt((t_vone *std_g.^2 + t_v * (reshape(curr_std_image, [size(centerview_video, 1), size(centerview_video, 2)])).^2)./(t_vone+t_v));
    save(sprintf('%s\\std_g%d_v%d.mat',std_savepath, g_id, cv_ind),'std_g','-v7.3');
    t_vone = t_vone + t_v;
end

save(sprintf('%s//shifts_wholeprocess.mat',outdir_p),'shifts','-v7.3');
t_ms = toc;
fprintf('measure shift process done and it takes %.2f secs\n',t_ms);
% save parameters
reg_param.reg_savepath = reg_savepath;
reg_param.std_savepath = std_savepath;
reg_param.bound = bound;
reg_param.option_r = option_r;
save(sprintf('%s\\param_reg.mat',outdir_p),'reg_param','-v7.3');
%% apply shift
load(sprintf('%s\\param_main.mat',outdir), 'main_param');
load(sprintf('%s\\param_debg.mat',outdir_p), 'debg_param');
load(sprintf('%s\\param_reg.mat',outdir_p), 'reg_param');
load(sprintf('%s\\shifts_wholeprocess.mat',outdir_p),'shifts');

debg_mul_video_savepath = debg_param.debg_mul_video_savepath;
reg_savepath = reg_param.reg_savepath;
std_savepath = reg_param.std_savepath;
savegroup = debg_param.savegroup;
maxIter = reg_param.maxIter;
rankdetrending = reg_param.rankdetrending;
option_r = reg_param.option_r;
bound = reg_param.bound;
view_array = view_config(main_param);
CoreNum_view = reg_param.CoreNum_view;

[~, ~, seq_ind] = Spiral_circle(main_param.Nnum,main_param.view_range);
seq_ind = seq_ind(end:-1:1);
t_regstd_start = clock; 
disp('--------------------------Apply shift for Nshift = 3--------------------------');

% apply shift for each angle

if isempty(gcp('nocreate'))
    p = parpool(CoreNum_view);
end
% for each view and parpool
parfor v_ind = 2 : length(view_array)
    tic;
    % apply shift
    v = seq_ind(v_ind);
    t_vone = 0;
    std_g = zeros(h, w);
    for g_id = 1: savegroup
        view_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff', debg_mul_video_savepath, v, g_id)));
        t_v = size(view_video,3);
        shifts_seg = shifts(t_vone + 1 : t_vone + t_v);
        view_video = apply_shifts(view_video, shifts_seg, option_r, bound/2, bound/2);
        
        saveastiff_n(v,im2uint16(view_video/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, v ,g_id));
        
        % calculate std
        
        if rankdetrending == 1
            [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(view_video, [], size(view_video,3)), maxIter);
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)), ...
                curr_bg_spatial(:), curr_bg_temporal);
        else
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)));
        end
        std_g = sqrt((t_vone * std_g.^2 + t_v * (reshape(curr_std_image, [size(view_video, 1), size(view_video, 2)])).^2)./(t_vone+t_v));
        parsave(sprintf('%s\\std_g%d_v%d.mat',std_savepath, g_id, v), std_g);
        t_vone = t_vone + t_v;
    end
    t_onereg = toc;
    fprintf('%d in %d view has been registered and std and it takes %.2f secs\n', v, length(view_array), t_onereg);
end
% delete(p);
t_regstd = etime(clock,t_regstd_start);
fprintf('resgistration process done and it takes %.2f secs\n',t_regstd); 
%% concate std
load(sprintf('%s\\param_main.mat',outdir),'main_param');
load(sprintf('%s\\param_debg.mat',outdir_p),'debg_param');
load(sprintf('%s\\param_reg.mat',outdir_p),'reg_param');
view_array = view_config(main_param);
std_savepath = reg_param.std_savepath;
first_file_name = main_param.first_file_name;
% std_concate
std_WDF = zeros(h, w, length(view_array));
for v = 1: length(view_array)
    load(sprintf('%s\\std_g%d_v%d.mat',std_savepath,debg_param.savegroup,v),'std_g');
    std_WDF(:,:,v) = std_g;
end
% save std
saveastiff(im2uint16(std_WDF / max(std_WDF(:))), sprintf('%s\\std_%s%d.tif', outdir_p, first_file_name,0));
%% Reconstruction with shift map
load(sprintf('%s\\param_main.mat',outdir),'main_param','realign_param','psf_param','recon_param');

first_file_name = main_param.first_file_name;
view_array = view_config(main_param);

recon_param.Nshift = realign_param.Nshift; %% Scanning Times = 3 x 3
recon_param.Nnum = realign_param.Nnum; %% 15 x 15 pixels behind each MicroLen
recon_param.upsampling = recon_param.Nnum / recon_param.Nshift / recon_param.ds;
recon_param.psf_layer_position = 1:1: recon_param.psf_end;

if psf_param.correction == 1
    recon_psfpath = sprintf('%s\\ROI%d%s',psf_param.reconpsf_perfix,patch_id, psf_param.reconpsf_surfix);
else
    recon_psfpath = psf_param.recon_psfpath;
end
frame = 0;
std_WDF = double(loadtiff(sprintf('%s\\std_%s%d.tif', outdir_p, first_file_name,frame)));


% if recon_param.shiftflag == 1
%     Nshift = realign_param.Nshift;
%     map_path = recon_param.map_path;
%     load(sprintf('%s\\map_iter5.mat',map_path),'map_wavshape','mapNshift');
%     map_wavshape = imresize(map_wavshape,[preprocess_param.wdf_h,preprocess_param.wdf_w]);
%     top_left = preprocess_param.margin_point(1,:);
%     map_wavshape = map_wavshape(tl(1)-top_left(1)+1 : rd(1)-top_left(1)+1, ...
%         tl(2)-top_left(2)+1 : rd(2)-top_left(2)+1,:,:);
%     
%     [coordinate1,coordinate2]=meshgrid(1:size(map_wavshape,2),1:size(map_wavshape,1));
%     std_WDF_shift = zeros(size(map_wavshape,1),size(map_wavshape,2),length(view_array));
%     for v = 1: length(view_array)
%         tmp = std_WDF(:,:,v);
%         tmp = interp2(coordinate1,coordinate2,tmp,coordinate1+map_wavshape(:,:,v,2),coordinate2+map_wavshape(:,:,v,1),'cubic',0);
%         std_WDF_shift(:,:,v) = tmp;
%     end
%     std_WDF = std_WDF_shift;
%     clear std_WDF_shift;
% end

% load psf
disp('---------------------------load psf---------------------------')


[psf, recon_param] = psf_load_module(recon_param,recon_psfpath,view_array);
std_WDF_up = imresize(std_WDF, ...
    [floor(size(std_WDF,1)*recon_param.upsampling/2)*2+1,...
    floor(size(std_WDF,2)*recon_param.upsampling/2)*2+1],'cubic');
[stdv_h, stdv_w, ~] = size(std_WDF_up); 


% reconstruct with normal deconvolution
disp('--------------------------reconstruction-----------------------');
% save file path
std_recon_savepath = sprintf('%s\\std_recon', outdir_p);
std_recon_name_perfix = sprintf('%s\\std_recon_volume', std_recon_savepath);
if ~exist(std_recon_savepath,'file')
    mkdir(std_recon_savepath);
end

std_volume = ones(stdv_h, stdv_w,size(psf,4));
std_volume = std_volume./sum(std_volume(:)).*sum(std_volume(:))./(size(std_volume,3)*size(std_volume,4));
std_volume = recon_module(recon_param, std_volume, psf, std_WDF_up, std_recon_name_perfix, view_array, frame);
std_volume_ori = std_volume;

recon_overlap = round(recon_param.upsampling.* main_param.overlap);
std_volume = std_volume(recon_overlap+1:end-recon_overlap,recon_overlap+1:end-recon_overlap,:);
[stdv_h, stdv_w, stdv_z] = size(std_volume);

% save std reconstruction result
imwriteTFSK(uint16(std_volume/max(std_volume(:))*65535),sprintf('%s\\final_recon.tiff',outdir_p));
save(sprintf('%s\\final_recon.mat',outdir_p),'std_volume','-v7.3');
save(sprintf('%s\\final_recon_ori.mat',outdir_p),'std_volume_ori','-v7.3');

% save parameters
recon_param.recon_overlap = recon_overlap;
recon_param.vsize = [stdv_h, stdv_w, stdv_z];
save(sprintf('%s\\param_recon.mat',outdir_p),'recon_param','-v7.3');


%% multi scale debg for video Nshift = 1
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'video_realign_param','psf_param','viddebg_param');
% multi scale debg param
realigndata_name_perfix = video_realign_param.realigndata_name_perfix;
viddebg_param.Nshift = video_realign_param.Nshift; 
viddebg_param.Nnum = video_realign_param.Nnum; 
viddebg_param.psf_layer_position = 1 : 4 : viddebg_param.psf_end;
view_array = view_config(main_param);
viddebg_param.view_array = view_array;
first_file_name = main_param.first_file_name;
maxframe = viddebg_param.maxframe;

% group in time lapse
savegroup = ceil(video_realign_param.valid_frame_num/maxframe);

% load psf
debg_psfpath = psf_param.debg_psfpath;
[depsf, viddebg_param] = psf_load_module(viddebg_param,debg_psfpath,view_array);
fprintf('load low resolution psf process done\n');

% multiscale detrending start
disp('---------------------------multiscale detrending Nshift = 1---------------------------');

% output file name
vid_recon_mul_savepath = sprintf('%s\\vid_recon_mul',outdir_p); 
vid_debg_mul_video_savepath = sprintf('%s\\vid_mul_video',outdir_p);

if ~exist(vid_recon_mul_savepath,'file') 
    mkdir(vid_recon_mul_savepath);
end
if ~exist(vid_debg_mul_video_savepath,'file') 
    mkdir(vid_debg_mul_video_savepath);
end


% bg ratio and z_range
view_video = zeros(h_ds, w_ds, length(view_array));
for v = 1: length(view_array)
    view_video(:,:,v) = double(imread(sprintf('%s_No%d.tif', realigndata_name_perfix, 0),...
        sub2ind([main_param.Nnum,main_param.Nnum],view_array{v}(2),view_array{v}(1)),...
        'PixelRegion',{[tl_ds(1),rd_ds(1)],[tl_ds(2),rd_ds(2)]}));
end
viddebg_param = bg_param_def(depsf,viddebg_param,view_video,vid_recon_mul_savepath,first_file_name);

%start de background
t_viddebg_start = clock;
[~, ~, seq_ind] = Spiral_circle(main_param.Nnum,main_param.view_range);
seq_ind = seq_ind(end:-1:1);

count = 0;
bg_ratio = gpuArray(single(viddebg_param.bg_ratio));

% debg frame by frame
for frame_i = 1 : video_realign_param.valid_frame_num
    view_video = zeros(h_ds, w_ds, length(view_array));
    % load image
    for v = 1: length(view_array)
        view_video(:,:,v) = double(imread(sprintf('%s_No%d.tif', realigndata_name_perfix, frame_i - 1),...
            sub2ind([main_param.Nnum,main_param.Nnum],view_array{v}(2),view_array{v}(1)),...
            'PixelRegion',{[tl_ds(1),rd_ds(1)],[tl_ds(2),rd_ds(2)]}));
    end
    % multiscale sub ground noise
    view_video = sub_mul_bg_module(depsf,viddebg_param,view_video,...
        vid_recon_mul_savepath,first_file_name,frame_i-1,bg_ratio);
    
    % save
    if mod(frame_i-1,maxframe) == 0
        count = count + 1;
        for v_ind = 1:length(view_array)
            v = seq_ind(v_ind);
            imwrite(im2uint16(view_video(:,:,v)/65535),strcat(vid_debg_mul_video_savepath,'\\',...
                'mul_video_view',num2str(v),'_g',num2str(count),'.tiff'));
        end
    else
        for v_ind = 1:length(view_array)
            v = seq_ind(v_ind);
            imwrite(im2uint16(view_video(:,:,v)/65535),strcat(vid_debg_mul_video_savepath,'\\',...
                'mul_video_view',num2str(v),'_g',num2str(count),'.tiff'),'WriteMode','append');
        end
    end
    tt = toc;
    fprintf('%d frame is done, take %.2fs\n',frame_i,tt);
end

t_viddebg = etime(clock,t_viddebg_start);
fprintf('Debg process done and it takes %.2f secs\n',t_viddebg);
clear bg_ratio;
clear depsf;

% save parameters
viddebg_param.savegroup = savegroup;
viddebg_param.vid_debg_mul_video_savepath = vid_debg_mul_video_savepath;
save(sprintf('%s\\param_viddebg.mat',outdir_p), 'viddebg_param', '-v7.3');

%% Registration for Nshift = 1
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'vidreg_param');
load(sprintf('%s\\param_viddebg.mat',outdir_p),'viddebg_param');
vid_debg_mul_video_savepath = viddebg_param.vid_debg_mul_video_savepath;
view_array = view_config(main_param);
disp('---------------------------Registration for Nshift = 1---------------------------');
% save path
vidreg_savepath = sprintf('%s\\vidreg_path', outdir_p);
if ~exist(vidreg_savepath,'file')
    mkdir(vidreg_savepath);
end

disp('--------------------------Measure shift for Nshift = 1--------------------------');
tic;
% measure shift
cv_ind = ceil(length(view_array)/2);
centerview_video = [];

for g_id = 1 : viddebg_param.savegroup
    load_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff',vid_debg_mul_video_savepath, cv_ind, g_id)));
    centerview_video = cat(3,centerview_video,load_video);
end

[d1,d2,~] = size(centerview_video);
[~, shifts, bound, option_r] = motion_correction(centerview_video, d1, d2, vidreg_param, outdir_p);

% apply shift to center view
t_vone = 0;
for g_id = 1:viddebg_param.savegroup
    view_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff', vid_debg_mul_video_savepath, cv_ind, g_id)));
    t_v = size(view_video,3);
    shifts_seg = shifts(t_vone + 1 : t_vone + t_v);
    view_video = apply_shifts(view_video, shifts_seg, option_r, bound/2, bound/2);
    
    saveastiff(im2uint16(view_video/65535),sprintf('%s\\vidreg_view_%d_g_%d.tiff', vidreg_savepath, cv_ind ,g_id));
    t_vone = t_vone + t_v;
end

save(sprintf('%s//shifts_vidprocess.mat',outdir_p),'shifts','-v7.3');
t_ms = toc;
fprintf('measure shift process done and it takes %.2f secs\n',t_ms);
% save parameters
vidreg_param.reg_savepath = vidreg_savepath;
vidreg_param.bound = bound;
vidreg_param.option_r = option_r;
save(sprintf('%s\\param_vidreg.mat',outdir_p),'vidreg_param','-v7.3');
%% apply shift for Nshift = 1
load(sprintf('%s\\param_main.mat',outdir), 'main_param');
load(sprintf('%s\\param_viddebg.mat',outdir_p), 'viddebg_param');
load(sprintf('%s\\param_vidreg.mat',outdir_p), 'vidreg_param');
load(sprintf('%s\\shifts_vidprocess.mat',outdir_p),'shifts');

vid_debg_mul_video_savepath = viddebg_param.vid_debg_mul_video_savepath;
savegroup = viddebg_param.savegroup;
vidreg_savepath = vidreg_param.reg_savepath;
option_r = vidreg_param.option_r;
bound = vidreg_param.bound;
view_array = view_config(main_param);
CoreNum_view = vidreg_param.CoreNum_view;

[~, ~, seq_ind] = Spiral_circle(main_param.Nnum,main_param.view_range);
seq_ind = seq_ind(end:-1:1);
t_vidreg_start = clock; 

disp('--------------------------Apply shift for Nshift = 1--------------------------');
% apply shift for each angle
if isempty(gcp('nocreate'))
    p = parpool(CoreNum_view);
end
% for each view and parpool
parfor v_ind = 2 : length(view_array)
    tic;
    % apply shift
    v = seq_ind(v_ind);
    view_video = [];
    t_vone = 0;
    for g_id = 1: savegroup
        view_video = double(loadtiff(sprintf('%s\\mul_video_view%d_g%d.tiff', vid_debg_mul_video_savepath, v, g_id)));
        t_v = size(view_video,3);
        shifts_seg = shifts(t_vone + 1 : t_vone + t_v);
        view_video = apply_shifts(view_video, shifts_seg, option_r, bound/2, bound/2);
        
        saveastiff_n(v,im2uint16(view_video/65535),sprintf('%s\\vidreg_view_%d_g_%d.tiff', vidreg_savepath, v ,g_id));
        t_vone = t_vone + t_v;
    end
    t_onereg = toc;
    fprintf('%d in %d view has been registered and it takes %.2f secs\n', v, length(view_array), t_onereg);
end
t_vidreg = etime(clock,t_vidreg_start);
fprintf('resgistration process done and it takes %.2f secs\n',t_vidreg); 

%% neuron segmentation
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'seed_param','realign_param','video_realign_param');
load(sprintf('%s\\param_recon.mat',outdir_p), 'recon_param');
% load mask
vessel_mask = loadtiff(sprintf('%s\\vessel_mask.tif',outdir_p));

seed_param.outdir = outdir_p;
seed_param.Nnum = main_param.Nnum;
seed_param.optical_psf_ratio = seed_param.per_slice_depth / seed_param.pixel_size;

seed_param.volume_threshold = [5 * seed_param.neuron_lateral_size.^2, 500* seed_param.neuron_lateral_size.^2];
seed_param.shell_radius = ceil(2 * seed_param.neuron_lateral_size);

seed_param.overlap = main_param.overlap;
seed_param.recon_overlap = recon_param.recon_overlap;
seed_param.margin = round(seed_param.neuron_lateral_size * 1.5);

seed_param.NN_ds = main_param.Nnum / realign_param.Nshift;
seed_param.vol_ds = recon_param.ds;
seed_param.wdf_ds = realign_param.Nshift/video_realign_param.Nshift;

recon_overlap = seed_param.recon_overlap;
overlap = seed_param.overlap;
margin = seed_param.margin;

% load std_volume
load(sprintf('%s\\final_recon.mat',outdir_p),'std_volume');
load(sprintf('%s\\final_recon_ori.mat',outdir_p),'std_volume_ori');

std_volume = double(std_volume);
std_volume_ori = double(std_volume_ori);

seed_param = zrange_def(std_volume,seed_param);
start_ind = 40; % can be manully changed
end_ind = 70;
P_th = 0.48;
seed_param.start_ind = start_ind;
seed_param.end_ind = end_ind;

max_value = 3050; % can be manully changed
min_value = 50;

seed_param.max_value = max_value;
seed_param.min_value = min_value;

std_volume = max(min(std_volume, max_value)-min_value, 0);
std_volume_ori = max(min(std_volume_ori, max_value)-min_value, 0);

if seed_param.mask3D == 0
    vessel_mask = vessel_mask(:,:,1);
else
    vessel_mask = vessel_mask(:,:,1 : 2: end);
end

vessel_mask(vessel_mask < P_th) = 0;
vessel_mask(vessel_mask >= P_th) = 1;

std_volume_cut = std_volume(:, :, start_ind : end_ind);
std_volume_cut = std_volume_cut / max(std_volume_cut(:));
std_volume_ori = std_volume_ori(:, :, start_ind : end_ind);
std_volume_ori = std_volume_ori / max(std_volume_ori(:));


% patch preparation
% --------------------patch preparation--------------------
[subpatch_info_array] = determine_patch_size_LFM(size(std_volume, 1), size(std_volume, 2), h, w, seed_param);


for global_patch_id = 1 : length(subpatch_info_array)% for lateral patches
    
    fprintf('patch %d\n', global_patch_id);
    curr_outdir = sprintf('%s\\subpatch_%d', outdir_p, global_patch_id);
    mkdir(curr_outdir);
    
    % volume preparation
    curr_patch_info = subpatch_info_array{global_patch_id};
    patch_volume = std_volume_cut(curr_patch_info.location(1, 1) : curr_patch_info.location(2, 1), ...
        curr_patch_info.location(1, 2) : curr_patch_info.location(2, 2), ...
        :);
    vessel_mask_cut = vessel_mask(curr_patch_info.location(1, 1) : curr_patch_info.location(2, 1), ...
        curr_patch_info.location(1, 2) : curr_patch_info.location(2, 2), ...
        :);
    patch_volume_ori = std_volume_ori(curr_patch_info.location(1, 1) + recon_overlap - margin: curr_patch_info.location(2, 1) + recon_overlap + margin, ...
        curr_patch_info.location(1, 2) + recon_overlap - margin: curr_patch_info.location(2, 2) + recon_overlap + margin, ...
        :);
    
    % neuron segmentation generation module
    center_array = [];
    disp('--------------------------Neuron segmentation--------------------------')
    % only keep central ones
    curr_seed_param = seed_param;
    curr_seed_param.outdir = curr_outdir;
    valid_seg = segmentation_module(patch_volume, patch_volume_ori, vessel_mask_cut, curr_seed_param);
    if ~isempty(valid_seg)
        for i = 1 : size(valid_seg, 1)
            center_array(i, :) = mean(valid_seg{i, 2}, 1);
        end
        % apply patches shifts
        reg_seg = valid_seg;
        reg_center = center_array;
        for k = 1 : size(valid_seg, 1)
            % note each seg has multiple small patches
            reg_seg{k, 2}(:, 1) = valid_seg{k, 2}(:, 1) + curr_patch_info.location(1, 1) - 1;
            reg_seg{k, 2}(:, 2) = valid_seg{k, 2}(:, 2) + curr_patch_info.location(1, 2) - 1;
            reg_seg{k, 2}(:, 3) = valid_seg{k, 2}(:, 3) + curr_seed_param.start_ind - 1;
            
            reg_center(:, 1) = center_array(:, 1) + curr_patch_info.location(1, 1) - 1;
            reg_center(:, 2) = center_array(:, 2) + curr_patch_info.location(1, 2) - 1;
            reg_center(:, 3) = center_array(:, 3) + curr_seed_param.start_ind - 1;
        end
        
        save(sprintf('%s\\neuron_seg_subpatch.mat', curr_outdir), 'reg_center', 'reg_seg', '-v7.3');
    else
        continue
    end
end
saveastiff(uint16(vessel_mask),sprintf('%s\\vessel_mask01.tiff',outdir_p));
save(sprintf('%s\\param_seed.mat',outdir_p),'seed_param','-v7.3');
%% CNMFE solve trace
load(sprintf('%s\\param_main.mat',outdir),'main_param', 'psf_param','realign_param', 'video_realign_param','demix_param');
load(sprintf('%s\\param_recon.mat',outdir_p), 'recon_param');
load(sprintf('%s\\param_viddebg.mat',outdir_p), 'viddebg_param');
load(sprintf('%s\\param_vidreg.mat',outdir_p), 'vidreg_param');
load(sprintf('%s\\param_seed.mat',outdir_p), 'seed_param');

load(sprintf('%s\\final_recon.mat',outdir_p),'std_volume');
view_array = view_config(main_param);
if psf_param.correction == 1
    recon_psfpath = sprintf('%s\\ROI%d%s',psf_param.reconpsf_perfix, patch_id, psf_param.reconpsf_surfix);
else
    recon_psfpath = psf_path.recon_psfpath;
end
[psf, ~] = psf_load_module(recon_param,recon_psfpath,view_array);

savegroup = viddebg_param.savegroup;
vidreg_savepath = vidreg_param.reg_savepath;
valid_frame_num = video_realign_param.valid_frame_num;
vsize = recon_param.vsize;

start_ind = seed_param.start_ind;
end_ind = seed_param.end_ind;

bg_iter = demix_param.bg_iter;
max_demixing_round = demix_param.max_demixing_round;
maxIter_NMF = demix_param.maxIter_NMF;
oasis_lambda = demix_param.oasis_lambda;
oasis_g = demix_param.oasis_g;
lambda_l0 = demix_param.lambda_l0;
frames_step = demix_param.frames_step;


% define wigners that is usefull
specify_wigner = zeros(5,1);
specify_wigner(1) = find(cellfun(@(x)all(x(:)==[ceil(main_param.Nnum / 2); ceil(main_param.Nnum / 2)]),view_array));
for i = 1 : 4
    [u, v] = ind2sub([2, 2], i);
    buf = [ceil(main_param.Nnum / 2) + (-1)^u * ceil(main_param.Nnum / 5),...
        ceil(main_param.Nnum / 2) + (-1)^v * ceil(main_param.Nnum / 5)];
    specify_wigner(i+1) = find(cellfun(@(x)all(x(:)==buf(:)),view_array));
end

[subpatch_info_array] = determine_patch_size_LFM(vsize(1), vsize(2), h, w, seed_param);
for global_patch_id = 1 : length(subpatch_info_array)% for lateral patches
    
    fprintf('patch %d\n', global_patch_id);
    curr_outdir = sprintf('%s\\subpatch_%d', outdir_p, global_patch_id);
    center_array = [];
    
    % volume preparation
    curr_patch_info = subpatch_info_array{global_patch_id};
    patch_volum_size_cut = [curr_patch_info.location(2, 1) - curr_patch_info.location(1, 1) + 1, ...
            curr_patch_info.location(2, 2) - curr_patch_info.location(1, 2) + 1, ...
            end_ind - start_ind + 1];
        
        
    % volume
    volume_img_cut = std_volume(curr_patch_info.location(1, 1):curr_patch_info.location(2, 1),...
        curr_patch_info.location(1, 2):curr_patch_info.location(2, 2), start_ind : end_ind);
    vol_ds = seed_param.vol_ds;
    wdf_ds = seed_param.wdf_ds;
    NN_ds = seed_param.NN_ds;
    ds_ratio = NN_ds * wdf_ds / vol_ds;
    volume_img_cut_ds = imresize(volume_img_cut,[curr_patch_info.wdf_ds_size(1),curr_patch_info.wdf_ds_size(2)]);
    psf_d = imresize(psf(:, :, :, start_ind:end_ind), [floor(size(psf,1)/ds_ratio/2)*2+1,floor(size(psf,2)/ds_ratio/2)*2+1],'cubic');
    wigner = zeros(size(volume_img_cut_ds,1),size(volume_img_cut_ds,2),size(specify_wigner, 1));

    for i = 1 : size(specify_wigner, 1)
        curr_vind = specify_wigner(i, 1);
        HXguess = prop_to_target_wigner(volume_img_cut_ds, psf_d, curr_vind);
        wigner(:,:,i) = HXguess;
    end
    
    frame_num = realign_param.num_rawdata;
    seg_stdvideo = zeros(size(wigner,1), size(wigner,2), size(specify_wigner, 1));
    group_num = ceil(frame_num/viddebg_param.maxframe);
    for j = 1 : size(specify_wigner, 1) % number of specified wigner
        curr_video = [];
        for g_id = 1 : group_num
            curr_video = cat(3,curr_video,loadtiffsub(sprintf('%s\\vidreg_view_%d_g_%d.tiff', vidreg_savepath, specify_wigner(j),g_id),...
                curr_patch_info.wdf_loc_ds(1,:),curr_patch_info.wdf_loc_ds(2,:)));
        end
        curr_video = curr_video(:,:,1:frame_num);
        if vidreg_param.rankdetrending == 1
            [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(curr_video, [], size(curr_video,3)), vidreg_param.maxIter);
            curr_std_image = compute_std_image(reshape(curr_video, [], size(curr_video,3)), ...
                curr_bg_spatial(:), curr_bg_temporal);
        else
            curr_std_image = compute_std_image(reshape(curr_video, [], size(curr_video,3)));
        end
        std_g = reshape(curr_std_image, [size(curr_video, 1), size(curr_video, 2)]);
        seg_stdvideo(:, :, j) = std_g;
        clear curr_video;
    end

    border_margin = 4;
    shift_block = zeros(2,size(specify_wigner, 1));
    [xx,yy] = meshgrid(1:size(wigner,2),1:size(wigner,1));
    for j = 1 : size(specify_wigner, 1)
        sub_wigner = wigner(border_margin + 1 : end - border_margin, border_margin + 1 : end - border_margin,j);
        sub_seg_stdvideo = squeeze(seg_stdvideo(:,:,j));
        corr_map = normxcorr2(sub_wigner,sub_seg_stdvideo);
        [shift_a,shift_b]=find(corr_map==max(corr_map(:)));
        shift_block(1,j)=shift_a(1)-size(sub_seg_stdvideo,1)+border_margin;
        shift_block(2,j)=shift_b(1)-size(sub_seg_stdvideo,2)+border_margin;
        if shift_block(1,j) > 3 || shift_block(1,j) < -3
            shift_block(1,j) = 2;
        end
        if shift_block(2,j) > 3 || shift_block(2,j) < -3
            shift_block(2,j) = 1;
        end
        wigner(:,:,j)=interp2(xx,yy,wigner(:,:,j),xx-shift_block(2,j),yy-shift_block(1,j),'cubic',0);
    end
    save(sprintf('%s\\shift_block.mat', curr_outdir),'shift_block','-v7.3');
    
    
    curr_seed_param = seed_param;
    curr_seed_param.outdir = curr_outdir;
    
    % load valid seg
    load(sprintf('%s\\final_filtering.mat', curr_outdir),'valid_seg_global_filt');
    valid_seg = valid_seg_global_filt;

    % if find(~cellfun(@isempty,valid_seg))
    if ~isempty(valid_seg)
        
        % seed generation
        disp('--------------------------Seed generation--------------------------')
        
        % define wigners that is usefull
        specify_wigner = zeros(5,1);
        specify_wigner(1) = find(cellfun(@(x)all(x(:)==[ceil(main_param.Nnum / 2); ceil(main_param.Nnum / 2)]),view_array));
        for i = 1 : 4
            [u, v] = ind2sub([2, 2], i);
            buf = [ceil(main_param.Nnum / 2) + (-1)^u * ceil(main_param.Nnum / 5),...
                ceil(main_param.Nnum / 2) + (-1)^v * ceil(main_param.Nnum / 5)];
            specify_wigner(i+1) = find(cellfun(@(x)all(x(:)==buf(:)),view_array));
        end
        
        % generate 2d initialized seed
        shell_radius = curr_seed_param.shell_radius;
        [S_init, S_shell_init,S_mask_init, S_shell_mask_init, valid_seg] = seed_generation_module(psf(:, :, :, start_ind:end_ind), valid_seg, ...
            patch_volum_size_cut, curr_patch_info, ...
            specify_wigner, shell_radius, curr_seed_param);
        
        % calcualte component center
        for i = 1 : size(valid_seg, 1)
            center_array(i, :) = mean(valid_seg{i, 2}, 1);
        end
        
        % prepare iteration video
        disp('--------------------------Patch video load--------------------------')
        [processed_video] = video_preparation_patch(savegroup, vidreg_savepath, valid_frame_num, specify_wigner, ...
            curr_patch_info, recon_param, curr_seed_param, S_init, shift_block);
        
        % background component initialization
        [bg_spatial_init, bg_temporal_init] = initialize_bg(processed_video, bg_iter);
        
        % temporal component initialization
        [T_init, T_shell_init] = initialize_T(processed_video, S_mask_init,  S_shell_mask_init);
        
        % change S shape
        % S
        S = cell_process_A(S_init, size(processed_video)); clear S_init
        S_mask_init= cell_process_A(S_mask_init, size(processed_video));
        S_mask_init = S_mask_init > 0;
        % S shell
        S_bg = cell_process_A(S_shell_init, size(processed_video)) ;  clear S_shell_init
        S_shell_mask_init = cell_process_A(S_shell_mask_init, size(processed_video));
        S_shell_mask_init = S_shell_mask_init > 0;
        
        save(sprintf('%s\\initialization.mat', curr_outdir), 'S', 'S_bg', 'S_mask_init', ...
            'S_shell_mask_init', 'T_init', 'T_shell_init', '-v7.3')
        
        % main iteration module
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
        
        % further core-shell demxing
        disp('--------------------------background subtraction--------------------------')
        tic
        neuron_trace_mat = zeros(size(T, 1), size(T, 2));
        deconvol_neuron_trace_mat = zeros(size(T, 1), size(T, 2));
        spike_mat = zeros(size(T, 1), size(T, 2));
        coefs_array = zeros(1, size(T, 1));
        
        % substract
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
        
        % apply patches shifts
        reg_seg = valid_seg;
        reg_center = center_array;
        for k = 1 : size(valid_seg, 1)
            % note each seg has multiple small patches
            reg_seg{k, 2}(:, 1) = valid_seg{k, 2}(:, 1) + curr_patch_info.location(1, 1) - 1;
            reg_seg{k, 2}(:, 2) = valid_seg{k, 2}(:, 2) + curr_patch_info.location(1, 2) - 1;
            reg_seg{k, 2}(:, 3) = valid_seg{k, 2}(:, 3) + curr_seed_param.start_ind - 1;
            
            reg_center(:, 1) = center_array(:, 1) + curr_patch_info.location(1, 1) - 1;
            reg_center(:, 2) = center_array(:, 2) + curr_patch_info.location(1, 2) - 1;
            reg_center(:, 3) = center_array(:, 3) + curr_seed_param.start_ind - 1;
        end
        
        save(fullfile(curr_outdir, ['registered.mat']),...
            'reg_center', 'reg_seg', 'neuron_trace_mat', 'deconvol_neuron_trace_mat','-v7.3');
    else
        continue
    end
    
end
%% stitch one patch
load(sprintf('%s\\param_recon.mat',outdir_p),'recon_param');
load(sprintf('%s\\param_seed.mat',outdir_p),'seed_param');
vsize = recon_param.vsize;
global_center = [];
global_seg = [];
global_trace = [];
global_trace_ori = [];
[subpatch_info_array] = determine_patch_size_LFM(vsize(1), vsize(2), h, w, seed_param);
for i = 1: length(subpatch_info_array)
    if exist(sprintf('%s\\subpatch_%d\\registered.mat',outdir_p,i),'file')
        load(sprintf('%s\\subpatch_%d\\registered.mat',outdir_p,i),'reg_center', 'reg_seg', 'neuron_trace_mat', 'deconvol_neuron_trace_mat');
        global_center = [global_center; reg_center];
        global_seg = [global_seg; reg_seg];
        global_trace = [global_trace; deconvol_neuron_trace_mat];
        global_trace_ori = [global_trace_ori; neuron_trace_mat];
    end
end
save(fullfile(outdir_p, ['global_output.mat']),...
    'global_center', 'global_seg', 'global_trace','global_trace_ori', '-v7.3');