function [valid_seg_global_filt] = segmentation_module(curr_volume, img, seed_param)
%% this module segment the reconstructed 3D image stack and output the neuron seed.
%  last update: 5/11/2021. YZ

%% Segmentation module
% segmentation parameters
% configure
% neuron_number = neuron_number ;
% seed_param.neuron_lateral_size = neuron_lateral_size;
% seed_param.local_constrast_th = local_constrast_th;
% seed_param.optical_psf_ratio = optical_psf_ratio;
% seed_param.overlap_ratio = overlap_ratio;
% seed_param.boundary_thickness = boundary_thickness;
% seed_param.discard_component_save = discard_component_save;
% seed_param.volume_threshold = [5 * seed_param.neuron_lateral_size.^2, 500* seed_param.neuron_lateral_size.^2];
outdir = seed_param.outdir;
mkdir(outdir)

ds_ratio = seed_param.NN_ds * seed_param.wdf_ds / seed_param.vol_ds;
% associate the reconstructions with central view.
estimate_patch_size = seed_param.estimate_patch_size;
pixel_size=seed_param.pixel_size;
per_slice_depth = seed_param.per_slice_depth;
% load current image    

% size
[size_h, size_w, size_z] = size(curr_volume);
psf_layer_position = 1 : size_z;
patch_info_array = determine_patch_size(size_h, size_w, estimate_patch_size);

%% loop for different patch
valid_seg_array = [];

for j = 1 : length(patch_info_array) % for lateral patches
    j
    curr_param = seed_param;
    curr_param.outdir = sprintf('%s\\seg_%d', curr_param.outdir, j);
    mkdir(curr_param.outdir)
    curr_patch_info = patch_info_array{j};

    % grab current patch
    patch_volume = curr_volume(curr_patch_info.location(1, 1) : curr_patch_info.location(2, 1), ...
                              curr_patch_info.location(1, 2) : curr_patch_info.location(2, 2), ...
                              :);
    % segmentation
    [valid_seg, discard_seg] = neuron_segmentation_module(patch_volume, curr_param);
    %%
    % apply the bias information
    if find(~cellfun(@isempty,valid_seg))
        for k = 1 : size(valid_seg, 1)
            % note each seg has multiple small patches
            valid_seg{k, 2}(:, 1) = valid_seg{k, 2}(:, 1) + curr_patch_info.location(1, 1);
            valid_seg{k, 2}(:, 2) = valid_seg{k, 2}(:, 2) + curr_patch_info.location(1, 2);
        end
        % concatenate the segmentation
        valid_seg_array = [valid_seg_array; valid_seg];
                  
    else   
        continue
    end
end
valid_seg_global = valid_seg_array;  

save(sprintf('%s\\valid_seg_global.mat', outdir), 'valid_seg_global')
%% boundary filter
% components in the boundary will be discarded
margin_dis =10;
exclude_ind = [];
for i = 1 : size(valid_seg_global, 1)
    curr_center = mean(valid_seg_global{i, 2}, 1);
    
    % boundary judge
    if abs(curr_center(1) - size_h) < margin_dis || abs(curr_center(2) - size_w) < margin_dis ||...
            abs(curr_center(1)) < margin_dis || abs(curr_center(2) ) < margin_dis 
        bound_flag = 0;
    else
        bound_flag = 1;
    end
      
    % calculate the result
    exclude_ind(i) = 1 - bound_flag;
end
valid_seg_global_bd_filt = valid_seg_global;
valid_seg_global_bd_filt(find(exclude_ind), :) = [];

save(sprintf('%s\\valid_seg_global_filt.mat', outdir), 'valid_seg_global_bd_filt')

% diagnose information: plot the global for 
plot_3D_distribution_mod_depth_irreg(valid_seg_global_bd_filt, [size_h, size_w], [1, 2] * pixel_size,...
                                                                psf_layer_position * per_slice_depth, 10, outdir,  'boundary_filter')

%% blood vessel mask
% run blood vessel extraction
% symmfilter = struct();
% symmfilter.sigma     = 8; % variance of DoG
% symmfilter.len       = 100; % rho, querying size
% symmfilter.sigma0    = 1; % blurress kernel
% symmfilter.alpha     = 0.2; % distance decreasing
% asymmfilter = false;
% 
% % down sample the image
% % image_d = imresize(image, 0.2);
% 
% % Apple BCOSFIRE filter
% % img = movie(:, :, 1);
% img = max(img(:)) - img;
% response_stack = BCOSFIRE_lfm(img, symmfilter, asymmfilter); % only 2d here
% figure, imshow(response_stack)
% 
% % calculate the mask
% se = strel('disk',1);
% response_stack_segm= imdilate(response_stack, se);
% response_stack_segm = response_stack_segm > 5e-2;
% %
% % response_stack_segm = imgaussfilt(response_stack_segm, 10);
% figure, imshow(response_stack_segm)
% saveastiff(im2uint16(response_stack_segm), sprintf('%s\\blood_vessel_mask.tiff', outdir));
% 
% response_stack_segm  = zeros(size(img));
response_stack_segm = img;
saveastiff(im2uint16(response_stack_segm), sprintf('%s\\blood_vessel_mask.tiff', outdir));
%% grab the center and spatial footprint in the central view stack

exclude_ind = 0;
for i = 1 : size(valid_seg_global_bd_filt, 1)
   i
   % load current segmentation
   curr_seg = valid_seg_global_bd_filt{i, 1};
   curr_center = valid_seg_global_bd_filt{i, 2};
   
   % create Ain 
   curr_A_in = generate_Ain(curr_seg, curr_center, size_h, size_w);
   
   % do the mapping
   curr_center_cut = mean(curr_center(:, 1 : 2), 1); % only xy
   curr_center_cut = curr_center_cut / ds_ratio;
   
    % concatenate

    % empty check
    if isempty(find(curr_A_in))
       exclude_ind(i) = 1; 
    end
    
    % vessel segment check
    if response_stack_segm(max(min(floor(curr_center_cut(1)), size(img,1)), 1), ...
            max(min(round(curr_center_cut(2)), size(img,2)), 1)) > 0 % in the mask
        exclude_ind(i) = 1; % set 0 to cancel other options
    end
    
end


% 3D segmentations
valid_seg_global_filt = valid_seg_global_bd_filt;
valid_seg_global_filt(find(exclude_ind), :) = [];
%
plot_3D_distribution_mod_depth_irreg(valid_seg_global_filt, [size_h, size_w], [1, 2] * pixel_size,...
                    psf_layer_position * per_slice_depth, 10, outdir, 'vessel_filter')
%
save(sprintf('%s\\final_filtering.mat', outdir),...
                                    'exclude_ind', 'valid_seg_global_filt')
end