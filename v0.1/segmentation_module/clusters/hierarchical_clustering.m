function [valid_seg, discard_seg] = hierarchical_clustering(segments, ...
    std_recon_img, boundary_thickness, optical_psf_ratio, overlap_ratio, volume_threshold, ...
    output_dir, discard_component_save)
%% in this file, I would use the hierarchical cluster (3D) of the generated
%  neuron segment, based on their 3D distribution.
%  use built-in command of matlab.
%  after clustering, do merging and discard the components that are too
%  large.

%  use text command to label the clusterred component
%  arrange different color
%  Input:
%       segments: per-plane segments
%       std_recon_img: raw segmented image
%       boundary_thickness: for plotting boundaries
%       output_dir: output dir
%       optical_psf_ratio: axial over lateral ratio

%  Output:
%       valid_seg: real valid segments (containes different depths)
%       discard_seg: real valid segments (containes different depths)

%  update: add boundary box and plot in the original image
%          add text to indicates the cluster number.
%          when do dilate, use boolean matrix for saving memory
%  to avoid large memory requirement, avoid the stack generation

%  last update: 3/20/2020. YZ
%% load the segments
% segmens = importdata('..\..\segmentation\raw_segments.mat');
% std_recon_img = loadtiff('..\..\example\with_deconv_recon_phasespace_theory_filter_LFM_sample_bg_size_450_450_51_N_neuron_100.tiff_iter_3.tiff');
% std_recon_img = double(std_recon_img);
% std_recon_img = std_recon_img / max(std_recon_img(:));
% std_recon_img = im2uint8(std_recon_img);
% size_h = 450;
% size_w = 450;
[size_h, size_w, size_z] = size(std_recon_img);
% boundary_thickness = 3; % thickness of boundary
%% generate patch center that distributed in 3D space
disp('cluster the segment')
% assign z-depth distance
% optical_psf_ratio = 1; % ratio of the optical psf size in axial/lateral

component_number = 0;
% merge different z
all_plane_seg = []; % 2 columns, 1st for segments, 2nd for center position
for i = 1 : size_z
    in_plane_seg = segments{i, 2};
    component_center = segments{i, 1};
    for j = 1 : size(component_center, 1)
        all_plane_seg{component_number + j, 1} = in_plane_seg{j};
        all_plane_seg{component_number + j, 2} = [component_center(j, :), i];
    end
    depth = (i - size_z / 2)  * optical_psf_ratio; % considering the enlarged psf
    component_center = [component_center, ones(size(component_center, 1), 1) * depth];
    position_matrx(component_number + 1: component_number + size(component_center, 1), 1 : 3) = component_center;
    component_number = component_number+ size(component_center, 1);
end

patch_size = size(all_plane_seg{1, 1}, 1);
%% do clustering
fprintf('\t\t hierarchical tree and cut... \n')

% calculate pairwise distance
% this distance can be changed, even based on the overlap.
D = pdist(position_matrx);

% build link, method as single
Z = linkage(D, 'single');

% cut based on the neuron size
T = cluster(Z, 'cutoff', overlap_ratio * patch_size, 'criterion', 'distance');

%% merging and discard.

% if the component is too large, discard it
N_seg = max(T);
valid_seg = [];
discard_seg = [];

img_volume = zeros(size_h, size_w, size_z);
% volume_threshold = patch_size^2 * 10 * 5;

valid_ind = 1;
discard_ind = 1;
try
textprogressbar('       valid segmentation merging ')
end
for i = 1 : N_seg
    try
    textprogressbar(i / N_seg * 100)
    end
    ind = find(T == i);
    
    segments_in_same_cluster = all_plane_seg(ind, 1);
    centers_in_same_cluster = cell2mat(all_plane_seg(ind, 2));
    % calculat the spatial span. now it has to introduce the square shape,
    % simple spot would not help.
    [dia_volume, volume_size] = volume_calculator(centers_in_same_cluster, ...
        patch_size, size_h, size_w, size_z);
    
    % judge the size
    if volume_size > volume_threshold(2) || volume_size < volume_threshold(1)
        % discard
        discard_seg{discard_ind, 1} = segments_in_same_cluster; % 2d segments
        discard_seg{discard_ind, 2} = centers_in_same_cluster; % 3d positions
        %         discard_seg{discard_ind, 3} = boolean(dia_volume); % volume, this might be too large for large dataset
        discard_ind = discard_ind + 1;
    else
        % record
        valid_seg{valid_ind, 1} = segments_in_same_cluster; % 2d segments
        valid_seg{valid_ind, 2} = centers_in_same_cluster; % 3d positions
        %         valid_seg{valid_ind, 3} = boolean(dia_volume); % volume, this might be too large for large dataset
        valid_ind = valid_ind + 1;
    end
end
try
textprogressbar(' done')
end

save(fullfile(output_dir,  'valid_segments.mat'), 'valid_seg', '-v7.3');
save(fullfile(output_dir,  'discard_segments.mat'), 'discard_seg', '-v7.3');
%% do fuse, use label function, XYCZT order
%  disable it for faster performance
%{1
img_valid = uint8(zeros(size_h, size_w, 3, size_z));
valid_num = size(valid_seg, 1);
colors = rand(valid_num, 3) * 0.5; % do a scale for clear background image

std_recon_img = double(std_recon_img);
std_recon_img = std_recon_img / max(std_recon_img(:));
std_recon_img = im2uint8(std_recon_img);

valid_in_std_recon_img = std_recon_img;
valid_in_std_recon_img = repmat(valid_in_std_recon_img, [1, 1, 1, 3]);
valid_in_std_recon_img = permute(valid_in_std_recon_img, [1, 2, 4, 3]);
% for bundary extraction
erode_se = strel('square', boundary_thickness * 2);
try
    textprogressbar('       create RGB label for valid segments')
end
for i = 1 : valid_num
    try
        textprogressbar(i / valid_num * 100)
    end
    %     volume = valid_seg{i, 3};
    [volume, ~] = volume_calculator( valid_seg{i, 2}, ...
        patch_size, size_h, size_w, size_z);
    for j = 1 : size_z
        
        % add segment label
        slice = squeeze(img_valid(:, :, :, j)); % 3 channel image
        bw_seg = logical(volume(:, :, j));
        slice = labeloverlay(slice, bw_seg, 'colormap', colors(i, :));
        
        % add component number
        ind = find(bw_seg);
        if ~isempty(ind)
            [sub_j, sub_i] = ind2sub([size_h, size_w], ind(1));
            try  % insertText() requires Computer Vision Toolbox
                slice = insertText(slice, [sub_i, sub_j], i,...
                    'Fontsize', 12, 'TextColor', colors(i, :) * 255, 'BoxOpacity', 0,...
                    'AnchorPoint', 'RightBottom'); % note this 255..
                % Too slow:
                %slice = AddTextToImage(slice, num2str(i), [sub_i, sub_j],...
                %            colors(i, :), 'Courier', 12);
            catch
            end
        end
        img_valid(:, :, :, j) = slice;
        % for debugging
        if j == 20
            %figure(101), imshow(slice), title(sprintf('component %d depth %d', i, j)), pause(0.01)
        end
        
        % add boundary box to the original image
        bw_seg_bound = bw_seg - imerode(bw_seg, erode_se);
        real_silce = squeeze(valid_in_std_recon_img(:, :, :, j));
        real_silce = labeloverlay(real_silce, bw_seg_bound, 'colormap', colors(i, :));
        
        % add component number
        if ~isempty(ind)
            try
                real_silce = insertText(real_silce, [sub_i, sub_j], i,...
                    'Fontsize', 12, 'TextColor', colors(i, :) * 255, 'BoxOpacity', 0,...
                    'AnchorPoint', 'RightBottom'); % note this 255..
                
                %slice = AddTextToImage(slice, num2str(i), [sub_i, sub_j],...
                %            colors(i, :), 'Courier', 12);
            catch
            end
        end
        valid_in_std_recon_img(:, :, :, j) = real_silce;
        if j == 20
            %figure(102), imshow(real_silce), title(sprintf('component %d depth %d', i, j)), pause(0.01)
        end
        
    end
end
try
    textprogressbar(' done')
end
hyperstack_write(fullfile(output_dir,'img_valid_segment.tif'), img_valid);
hyperstack_write(fullfile(output_dir, 'valid_in_std_recon_img.tif'), valid_in_std_recon_img);
%}
%% Plot discarded components: do fusion, use label function, XYCZT order
if discard_component_save
    img_discard = uint8(zeros(size_h, size_w, 3, size_z));
    discard_num = size(discard_seg, 1);
    colors = rand(discard_num, 3);
    
    discard_in_std_recon_img = std_recon_img;
    discard_in_std_recon_img = repmat(discard_in_std_recon_img, [1, 1, 1, 3]);
    discard_in_std_recon_img = permute(discard_in_std_recon_img, [1, 2, 4, 3]);
    % for bundary extraction
    erode_se = strel('square', boundary_thickness * 2);
    textprogressbar('       create RGB label for discard segments')
    for i = 1 : discard_num
        textprogressbar(i / discard_num * 100)
        [volume, ~] = volume_calculator( discard_seg{i, 2}, ...
            patch_size, size_h, size_w, size_z);
        for j = 1 : size_z
            % colorful area label
            slice = squeeze(img_discard(:, :, :, j)); % 3 channel image
            bw_seg = volume(:, :, j);
            slice= labeloverlay(slice, bw_seg, 'colormap', colors(i, :));
            
            % add component number
            % add component number
            ind = find(bw_seg);
            if ~isempty(ind)
                [sub_j, sub_i] = ind2sub([size_h, size_w], ind(1));
                try
                    slice = insertText(slice, [sub_i, sub_j], i,...
                        'Fontsize', 12, 'TextColor', colors(i, :) * 255, 'BoxOpacity', 0,...
                        'AnchorPoint', 'RightBottom'); % note this 255..
                catch
                end
            end
            
            img_discard(:, :, :, j) = slice;
            % for debugging
            if j == 20
                figure(101), imshow(slice), title(sprintf('component %d depth %d', i, j)), pause(0.01)
            end
            
            % add boundary box to the original image
            bw_seg_bound = bw_seg - imerode(bw_seg, erode_se);
            real_silce = squeeze(discard_in_std_recon_img(:, :, :, j));
            real_silce = labeloverlay(real_silce, bw_seg_bound, 'colormap', colors(i, :));
            
            % add component number
            if ~isempty(ind)
                try
                    real_silce = insertText(real_silce, [sub_i, sub_j], i,...
                        'Fontsize', 12, 'TextColor', colors(i, :) * 255, 'BoxOpacity', 0,...
                        'AnchorPoint', 'RightBottom'); % note this 255..
                catch
                end
            end
            discard_in_std_recon_img(:, :, :, j) = real_silce;
            if j == 20
                %figure(102), imshow(real_silce), title(sprintf('component %d depth %d', i, j)), pause(0.01)
            end
        end
    end
    textprogressbar('done')
    % write hyperstack
    hyperstack_write(fullfile(output_dir, [datestr(now, 'YYmmddTHHMM') '_img_discard_segment.tif']), img_discard);
    hyperstack_write(fullfile(output_dir, [datestr(now, 'YYmmddTHHMM') '_discard_in_std_recon_img.tif']), discard_in_std_recon_img);
    
else
    discard_seg = [];
end
end
