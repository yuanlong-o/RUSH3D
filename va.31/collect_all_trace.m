%% This file used to collect all trace and concate to one std image
% Author. MW
clc;clear;
close all;
addpath(genpath('utility'));
%% ROI 
ROI = 12;
bx = 3;
by = 4;

% file path
filepath = 'E:\20230314_64';
load(sprintf('%s\\patch_info.mat',filepath),'patch_info_array');
load(sprintf('%s\\param_main.mat',filepath),'main_param', 'realign_param', 'video_realign_param', 'recon_param','seed_param');

% important parameter will be used
final_param.pixel_size = seed_param.pixel_size;
final_param.per_slice_depth = seed_param.per_slice_depth;
final_param.neuron_lateral_size = seed_param.neuron_lateral_size;
final_param.local_constrast_th = seed_param.local_constrast_th;
final_param.overlap = main_param.overlap;
final_param.margin_point = main_param.margin_point;
final_param.Nshift1 = realign_param.Nshift;
final_param.Nshift2 = video_realign_param.Nshift;
final_param.recon_ds = recon_param.ds;


% trace collect
whole_center = [];
whole_seg = [];
whole_trace = [];
whole_trace_ori = [];
std_patch_info_array = [];
tmp_patch_loc = [1,1];
for patch_id = 1 : ROI
    fprintf('%d patch loading ...\n',patch_id);
    load(sprintf('%s\\patch_%d\\global_output.mat',filepath,patch_id),'global_center', 'global_seg', 'global_trace','global_trace_ori');
    load(sprintf('%s\\patch_%d\\param_recon.mat',filepath,patch_id),'recon_param');
    h = recon_param.vsize(1);
    w = recon_param.vsize(2);
    top_left = tmp_patch_loc;
    right_down = tmp_patch_loc + [h-1, w-1];
    std_patch_info.location = [top_left; right_down]; 
    std_patch_info.size = [h, w];
    std_patch_info_array{patch_id} = std_patch_info;
    tmp_seg = global_seg;
    tmp_center = global_center;
    for k = 1 : size(global_seg, 1)
        % note each seg has multiple small patches
        tmp_seg{k, 2}(:, 1) = global_seg{k, 2}(:, 1) + tmp_patch_loc(1) - 1;
        tmp_seg{k, 2}(:, 2) = global_seg{k, 2}(:, 2) + tmp_patch_loc(2) - 1;
        tmp_seg{k, 2}(:, 3) = global_seg{k, 2}(:, 3);
        
        tmp_center(:, 1) = global_center(:, 1) + tmp_patch_loc(1) - 1;
        tmp_center(:, 2) = global_center(:, 2) + tmp_patch_loc(2) - 1;
        tmp_center(:, 3) = global_center(:, 3);
    end
    
    whole_seg = [whole_seg;tmp_seg];
    whole_center = [whole_center;tmp_center];
    whole_trace = [whole_trace;global_trace];
    whole_trace_ori = [whole_trace_ori;global_trace_ori];
    
    if mod(patch_id,by) == 0
        tmp_patch_loc = [tmp_patch_loc(1)+h,1];
    else
        tmp_patch_loc = tmp_patch_loc + [0, w];
    end
   
end

save(sprintf('%s\\wholebrain_output.mat',filepath),...
    'whole_center', 'whole_seg', 'whole_trace','whole_trace_ori','final_param','-v7.3');
%% image tiff stitch
save(sprintf('%s\\std_patch_info.mat',filepath),'std_patch_info_array', '-v7.3');
whole_brain_3d = zeros(std_patch_info_array{ROI}.location(2,1),std_patch_info_array{ROI}.location(2,2),recon_param.vsize(3));
for patch_id = 1:ROI
    loc = std_patch_info_array{patch_id}.location;
    whole_brain_3d(loc(1,1):loc(2,1),loc(1,2):loc(2,2),:) = loadtiff(sprintf('%s\\patch_%d\\final_recon.tiff',filepath,patch_id));
end
saveastiff(uint16(whole_brain_3d),sprintf('%s\\whole_brain_3d.tiff',filepath));
%% plot neuron location
scatter(whole_center(:,2),whole_center(:,1),2,'filled');
axis([1,size(whole_brain_3d,2),1,size(whole_brain_3d,1)]);
axis off;
hold off;
axis tight;