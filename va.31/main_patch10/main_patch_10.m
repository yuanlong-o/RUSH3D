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
patch_id = 10; %1-12
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

XIAO