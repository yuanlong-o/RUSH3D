
function raw_segments = morphology_segment(recon_stack, neuron_number, neuron_lateral_size,...
    local_constrast_th, output_dir)
%% 3D segmentation. with two-parts
%  part 1: 3d morphology segmentation
%  Input:
%       recon_stack: input reconstructed stack
%       neuron_number: estimate neuron number
%       neuron_lateral_size: lateral neuron size
%       neuron_axial_size: axial neuron size
%       local_constrast_th: estimate neuron contrast compared to
%                        enviornment
%       output_dir: output directory
%       plot_flag: for debug
%  Output:
%       raw_segments: plane-by-plane segmentation
%  last update: 2/13/2020. YZ

%% load example 3D data
% recon_stack = loadtiff('example\with_deconv_recon_phasespace_theory_filter_LFM_sample_bg_size_450_450_51_N_neuron_100.tiff_iter_3.tiff');
% recon_stack = double(recon_stack);

[size_h, size_w, size_z] = size(recon_stack);
%% parameter
% neuron_number = 100;
% neuron_lateral_size = 10; % in pixel, diameter
% neuron_axial_size = 10; % in pixel, diameter

% 2d neuron template
% generate neuron template
gSig = [neuron_lateral_size / 2, neuron_lateral_size / 2]; % sigma
for i = 1 : 2
    gSiz_neuron(i) = ceil(3 * gSig(i) + 1);
    if mod(gSiz_neuron(i), 2) == 0
        gSiz_neuron(i) = gSiz_neuron(i) + 1;
    end
end
x = (1 : gSiz_neuron(1)) - ceil(gSiz_neuron(1) / 2);
y = (1 : gSiz_neuron(2)) - ceil(gSiz_neuron(2) / 2);
[Y, X] = meshgrid(y, x);
neuron_temp = exp(-((X / gSig(1)).^2 + (Y / gSig(2)).^2));
neuron_temp = neuron_temp / sum(neuron_temp(:));
% define local area template
box = round(gSiz_neuron(1) * 1.5);
box_half = ceil(box / 2);
box = box_half * 2  - 1;

local_temp = ones(box, box);
local_temp(box_half - floor(gSiz_neuron(1) / 2) : box_half + floor(gSiz_neuron(1) / 2), ...
    box_half - floor(gSiz_neuron(2) / 2) : box_half +floor(gSiz_neuron(2) / 2)) = 0;
local_temp = local_temp / sum(local_temp(:));

figure, imagesc(neuron_temp)
figure, imagesc(local_temp)
%% part 1: finding neurons frame-by-frame 
raw_segments = cell(size_z, 2);
simu_segment_volume = zeros(size_h, size_w, size_z);
real_segment_volume = zeros(size_h, size_w, size_z);
residue_voume = zeros(size_h, size_w, size_z);


% local_constrast_th = 1.1;
total_find = 0;
plot_flag = false;

textprogressbar('per depth 2D segmentation')
for i = 1 : size_z
    textprogressbar(i / size_z * 100)
    slice = recon_stack(:, :, i);
    slice = slice - min(slice(:));
    [simulate_segment, real_segment, record_matrix, res, update_K, component_lib] = ...
    suit_match_segment_2D(slice, neuron_temp, local_temp,gSiz_neuron(1), ...
    neuron_number,  local_constrast_th, i, plot_flag);

    raw_segments{i, 1} = record_matrix; % center of segments
    raw_segments{i, 2} = component_lib; % recorded segments
    
    simu_segment_volume(:, :, i) = simulate_segment;
    residue_voume(:, :, i) = res;
    real_segment_volume(:, :, i) = real_segment;
    total_find = total_find + update_K;
end
textprogressbar('done')
%% 
saveastiff(im2uint16(simu_segment_volume/ max(simu_segment_volume(:))), ...
    sprintf('%s/simu_segment_volume.tiff', output_dir))
saveastiff(im2uint16(residue_voume / max(residue_voume(:))), ...
    sprintf('%s/residue_voume.tiff', output_dir))
saveastiff(im2uint16(real_segment_volume / max(real_segment_volume(:))), ...
    sprintf('%s/real_segment_volume.tiff', output_dir))

save(sprintf('%s/raw_segments.mat', output_dir), 'raw_segments')
