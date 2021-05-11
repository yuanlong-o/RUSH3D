function [valid_seg, discard_seg] = neuron_segmentation_module(volume, SI)
%% module for 3D neuron reconstruction.
%   Input
%   volume: reconstructed 3D stack
%   SI: configuration structure
%       SI.neuron_number: maximum per depth neuron number
%       SI.neuron_lateral_size: lateral neuron size
%       SI.local_constrast_th: neuron/background contrast
%       SI.optical_psf_ratio : psf extend ratio
%       SI.overlap_ratio: merging ratio based on patch size
%       SI.boundary_thickness: for plot; show a box for each neuron
%       SI.discard_component_save: flag, if one will save the discard
%       component
%       SI.outdir: output directory

%   Output
%   valid_seg: cell array, record each valid component
%   discard_seg: cell array, record discarded component (as candicates of neurons)

%   last update: 9/8/2020. YZ

%% parser
neuron_number = SI.neuron_number;% maximum per depth neuron number
neuron_lateral_size = SI.neuron_lateral_size;
local_constrast_th = SI.local_constrast_th;

optical_psf_ratio = SI.optical_psf_ratio;
overlap_ratio = SI.overlap_ratio;
boundary_thickness = SI.boundary_thickness; % for plotting
discard_component_save = SI.discard_component_save;

volume_threshold = SI.volume_threshold;
outdir = SI.outdir;

d1 = size(volume, 1);
d2 = size(volume, 2);
d3 = size(volume, 3);
%% Segmentation
% neuron_number = 200;
% neuron_lateral_size = 15;
% local_constrast_th = 1.3;

%% call segmentation function
options.false_threshold = 100;
options.local_constrast_th = local_constrast_th;
options.d1 = d1;
options.d2 = d2;
options.min_threshold = 0.5;
options.gSig = neuron_lateral_size;    % width of the gaussian kernel approximating one neuron
options.gSiz = round(neuron_lateral_size * 2.5);    % average size of filter size
options.ssub = 1;   

options.min_v_search = 0.05;
options.seed_method  = 'other';
options.min_pixel = neuron_lateral_size^2 * 0.8;
options.center_psf = [];
options.spatial_constraints.connected =  true;


%% preparing
volume = volume / max(volume(:));
% volume = imadjustn(volume);

Cn_stack = zeros(d1, d2, d3);
raw_segments = cell(d3, 2);
simu_segment_volume = zeros(d1, d2, d3);
residue_voume = zeros(d1, d2, d3);



for i = 1 : d3
    i
    % input
    slice = volume(:, :, i);
    [results, center, Cn] = greedyROI_endoscope_summary_shape_aware(slice(:), ...
                                            neuron_number, options, false);
                                        
    Cn_stack(:, :, i) = Cn;
    center_stack{i} = center;
    
    % record raw_segments
    for j = 1 : size(center, 1)
        buf = results.Ain(:, j);
        buf = reshape(full(buf), d1, d2);
        patch = buf(max(center(j, 1) - options.gSig, 1) : min(center(j, 1) + options.gSig, d1), ...
            max(center(j, 2) - options.gSig, 1) : min(center(j, 2) + options.gSig, d2));
        component_lib{j} = patch;
    end
    
    if  isempty(center)
        raw_segments{i, 1} =[];
        raw_segments{i, 2} = [];
        
        simu_segment_volume(:, :, i) = zeros(d1, d2);
        residue_voume(:, :, i) = slice;
    else
        raw_segments{i, 1} = center; % center of segments
        raw_segments{i, 2} = component_lib; % recorded segments
        simulate_segment = reshape(full(sum(results.Ain, 2)), d1, d2);
    
        simu_segment_volume(:, :, i) = simulate_segment;
        residue_voume(:, :, i) = slice - simulate_segment;
    end
    % record dummy volume

    close all
end          
saveastiff(im2uint16(simu_segment_volume/ max(simu_segment_volume(:))), ...
    sprintf('%s/simu_segment_volume.tiff', outdir ))
saveastiff(im2uint16(residue_voume / max(residue_voume(:))), ...
    sprintf('%s/residue_voume.tiff', outdir ))
saveastiff(im2uint16(Cn_stack / max(Cn_stack(:))), ...
    sprintf('%s/Cn_stack.tiff', outdir ))
save(sprintf('%s/raw_segments.mat', outdir), 'raw_segments')
%% Clustering of segments
% optical_psf_ratio = 3;
% overlap_ratio = 0.5;
% boundary_thickness = 3; % for plotting
% discard_component_save = false;
% volume_threshold = [(neuron_lateral_size * 1.5)^2 * 5, ... % maximum/minimum neuron size
%     (neuron_lateral_size * 1.5)^2 * 500];
if  numel(find(~cellfun(@isempty, raw_segments))) > 4
      [valid_seg, discard_seg] = hierarchical_clustering(raw_segments, ...
        volume, boundary_thickness, optical_psf_ratio, overlap_ratio, ...
        volume_threshold, SI.outdir, discard_component_save); 
else
   
	valid_seg = cell(1, 2);
    discard_seg = cell(1, 2);
end

