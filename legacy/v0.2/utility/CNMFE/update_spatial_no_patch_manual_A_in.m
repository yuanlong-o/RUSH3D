function [obj, A_maunal_update, empty_ind ]= update_spatial_no_patch_manual_A_in(Y_mat, A_maunal, obj, update_sn)
%% update the the spatial components for all neurons
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.
%   update_sn:  boolean, update noise level for each pixel

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%  modified by YZ. add a source detection for more robust performance.
%  last update: 10/15/2020. YZ
%% process parameters

try
    % map data
    d1 = size(Y_mat, 1);
    d2 = size(Y_mat, 2);
    T = size(Y_mat, 3);
    obj.options.d1 = d1;
    obj.options.d2 = d2;
    
	patch_pos = [[1; d1], [1; d2]];    % patch position
    block_pos = [[1; d1], [1; d2]];    % patch position
catch
    error('No data file selected');
end
fprintf('\n-----------------UPDATE SPATIAL---------------------------\n');
% frames to be loaded
frame_range = obj.frame_range;
T = diff(frame_range) + 1;

% threshold for detecting large residuals
% thresh_outlier = obj.options.thresh_outlier;
% update sn or not
if ~exist('update_sn', 'var')||isempty(update_sn)
    update_sn = false; %don't save initialization procedure
end
% options
options = obj.options;
bg_model = options.background_model;
bg_ssub = options.bg_ssub;
method = options.spatial_algorithm;

%% determine search location for each neuron
search_method = options.search_method;
if strcmpi(search_method, 'dilate')
    obj.options.se = [];
end
% pre-remove all zeros component

i = 1;
num_neuron = size(obj.A, 2);
IND = A_maunal>0; % not A_manual is used for forcing a searching area
A_maunal_update = A_maunal;
empty_ind = zeros(1, num_neuron);
while i <= num_neuron
    num_neuron = size(obj.A, 2);
%     i
    if isempty(find(obj.A(:, i), 1))
%        obj.delete(i) 
        obj = CNMFE_delete(obj, i);
        IND(:, i) = []; % remember also delete IND
        A_maunal_update(:, i) = [];
        empty_ind(i) = 1;
    else
        i = i + 1;
    end
    
end
% IND = sparse(logical(determine_search_location(obj.A, search_method, options)));
% IND = true(d1 * d2, num_neuron); % force search

%% identify existing neurons within each patch

tmp_patch = patch_pos;
tmp_block = block_pos;
% find the neurons that are within the block
mask = zeros(d1, d2);
mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
mask(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)) = 2;
% find neurons within the patch (all neurons)
ind = find(reshape(mask(:)==2, 1, [])* full(double(IND))>0); % find here is time consuming, since it requires full
A= obj.A((mask>0), ind);
IND_search = IND(mask==2, ind);
sn = obj.P.sn(mask==2);
C = obj.C(ind, :);
ind_neurons = ind;    % indices of the neurons within each patch
with_neuron = ~isempty(ind);

if strcmpi(bg_model, 'ring')
    ind = find(reshape(mask(:)==1, 1, [])* full(obj.A_prev)>0); % find here is time consuming
    A_prev= obj.A_prev((mask>0), ind);
    C_prev = obj.C_prev(ind, :);
end

if update_sn
    sn_new = sn;
end
%% prepare for the variables for computing the background.
bg_model = obj.options.background_model;
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;

%% start updating spatial components

% prepare for updating model variables
tmp_patch = patch_pos;     %[r0, r1, c0, c1], patch location
if strcmpi(bg_model, 'ring')
    tmp_block = block_pos;
else
    tmp_block = patch_pos;
end
A_patch = A;
C_patch = C;                % previous estimation of neural activity
IND_patch = IND_search;
nr = diff(tmp_patch(1:2)) + 1;
nc = diff(tmp_patch(3:4)) + 1;

% use ind_patch to indicate pixels within the patch
ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;

% get data
if strcmpi(bg_model, 'ring')
    % including areas outside of the patch for recorving background
    % in the ring model
    Ypatch = Y_mat(:, :, frame_range(1) : frame_range(2));;
else
    error('only ring mode supported')
end
[nr_block, nc_block, ~] = size(Ypatch);

% get the noise level
sn_patch = sn;

% get background
if strcmpi(bg_model, 'ring')
    A_patch_prev = A_prev;
    C_patch_prev = C_prev;
    W_ring = W;
    b0_ring = b0;
    Ypatch = reshape(Ypatch, [], T);
    tmp_Y = double(Ypatch)-A_patch_prev*C_patch_prev;

    if bg_ssub==1
        Ypatch = bsxfun(@minus, double(Ypatch(ind_patch,:))- W_ring*tmp_Y, b0_ring-W_ring*mean(tmp_Y, 2));
    else
        % get the dimension of the downsampled data
        [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
        % downsample data and reconstruct B^f
        temp = reshape(bsxfun(@minus, tmp_Y, mean(tmp_Y, 2)), nr_block, nc_block, []);
        temp = imresize(temp, 1./bg_ssub);
        Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
        Bf = imresize(Bf, [nr_block, nc_block]);
        Bf = reshape(Bf, [], T);

        Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring);
    end
elseif strcmpi(bg_model, 'nmf')
    b_nmf = b;
    f_nmf = f;
    Ypatch = double(reshape(Ypatch, [], T))- b_nmf*f_nmf;
else
    b_svd = b;
    f_svd = f;
    b0_svd = b0;
    Ypatch = double(reshape(Ypatch, [], T)) - bsxfun(@plus, b_svd*f_svd, b0_svd);
end

% using HALS to update spatial components
if update_sn
    sn_patch = GetSn(Ypatch);
    sn_new = reshape(sn_patch, nr, nc);
end

A_patch = A_patch(ind_patch, :);

%         temp = HALS_spatial(Ypatch, A_patch, C_patch, IND_patch, 3);
if strcmpi(method, 'hals')
    temp = HALS_spatial(Ypatch, A_patch, C_patch, IND_patch, 3);
elseif strcmpi(method, 'hals_thresh')
    temp = HALS_spatial_thresh(Ypatch, A_patch, C_patch, IND_patch, 3, sn_patch);
elseif strcmpi(method, 'lars')
    temp = lars_spatial(Ypatch, A_patch, C_patch, IND_patch, sn_patch);
    %         elseif strcmpi(method, 'nnls_thresh')&&(~isempty(IND_patch))
    %             temp = nnls_spatial_thresh(Ypatch, A_patch, C_patch, IND_patch, 5, sn_patch);
else
    temp = nnls_spatial(Ypatch, A_patch, C_patch, IND_patch, 20);
end
%         A_new = tmp_obj.post_process_spatial(reshape(full(temp), nr, nc, []));
A_new = full(temp);


%% collect results
K = size(obj.A, 2);
A_ = zeros(d1, d2, K);
fprintf('Collect results from all small patches...\n');

A_patch = A_new;
tmp_pos = patch_pos;
nr = diff(tmp_pos(1:2))+1;
nc = diff(tmp_pos(3:4))+1;
ind_patch = ind_neurons;
for m=1:length(ind_patch)
    k = ind_patch(m);
    A_(tmp_pos(1):tmp_pos(2), tmp_pos(3):tmp_pos(4), k) = reshape(A_patch(:, m), nr, nc, 1);
end

A_new = sparse( CNMFE_reshape(obj,A_, 1)); % still sparse. but why?
if update_sn
    obj.P.sn = sn_new;
end
%% post-process results
fprintf('Post-process spatial components of all neurons...\n');
obj.A = post_process_spatial(obj, CNMFE_reshape(obj, A_new, 2)); % note this reshape will use full!
% obj.A = A_new;
fprintf('Done!\n');
fprintf('Done!\n');

%% upadte b0
if strcmpi(bg_model, 'ring')
    fprintf('Update the constant baselines for all pixels..\n');
    obj.b0_new = obj.P.Ymean-CNMFE_reshape(obj, obj.A*mean(obj.C,2), 2); %-obj.reconstruct();
    fprintf('Done!\n');
end

