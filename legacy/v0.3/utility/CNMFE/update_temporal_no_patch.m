function obj = update_temporal_no_patch(Y_mat, obj, use_c_hat)
%% update the the temporal components for all neurons
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.
%   use_c_hat: use the previous estimation of C

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

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
fprintf('\n-----------------UPDATE TEMPORAL---------------------------\n');

% frames to be loaded
frame_range = obj.frame_range;
T = diff(frame_range) + 1;


% use previous estimation of C
if ~exist('use_c_hat', 'var') || isempty(use_c_hat)
    use_c_hat = true;
end
% options
options = obj.options;
bg_model = options.background_model;
maxIter = options.maxIter;

%% identify existing neurons within each patch
if strcmpi(bg_model, 'ring')
    tmp_block = block_pos;
else
    tmp_block = patch_pos;
end
% find the neurons that are within the block
mask = zeros(d1, d2);
mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
ind = (reshape(mask(:), 1, [])* obj.A>0);
A= obj.A(logical(mask), ind);
sn = obj.P.sn(logical(mask));
C = obj.C(ind, :);
ind_neurons = find(ind);    % indices of the neurons within each patch

if strcmpi(bg_model, 'ring')
    ind = find(reshape(mask(:)==1, 1, [])* full(obj.A_prev)>0);
    A_prev= obj.A_prev((mask>0), ind);
    C_prev = obj.C_prev(ind, :);
end

%% prepare for the variables for computing the background.
bg_model = obj.options.background_model;
bg_ssub = obj.options.bg_ssub;
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;

%% start updating temporal components
deconv_flag = obj.options.deconv_flag;
if options.deconv_flag
    deconv_options = obj.options.deconv_options;
else
    deconv_options = [];
end

% no neurons within the patch
tmp_patch = patch_pos;     %[r0, r1, c0, c1], patch location
if strcmpi(bg_model, 'ring')
    tmp_block = block_pos;
else
    tmp_block = patch_pos;
end

C_patch = C;                % previous estimation of neural activity
A_patch = A;

% use ind_patch to indicate pixels within the patch
ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;

% get data
if strcmpi(bg_model, 'ring')
    % including areas outside of the patch for recorving background
    % in the ring model
    Ypatch = Y_mat(:, :, frame_range(1) : frame_range(2));
else
    error('only ring mode supported')
end
[nr_block, nc_block, ~] = size(Ypatch);

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

        Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring); % minus bg
    end
else
    error('only ring mode supported')
end

% using HALS to update temporal components
if ~use_c_hat
    [AA, C_raw_new] = fast_temporal(Ypatch, A_patch(ind_patch,:));
else
    [~, C_raw_new] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, maxIter, deconv_options); % where is bg?
    AA= sum(A_patch(ind_patch,:).^2, 1);
end

%% collect results
K = size(obj.C, 1);
C_new = zeros(K, T);
aa = zeros(K, 1);
fprintf('Collect results from all small patches...\n');

C_raw_patch = C_raw_new;
ind_patch = ind_neurons;
aa_patch = AA;
for m=1:length(ind_patch)
    k = ind_patch(m);
    C_new(k, :) = C_new(k, :) + C_raw_patch(m,:) * aa_patch(m);
    aa(k) = aa(k) + aa_patch(m);
end

aa(aa==0) = 1;
obj.C_raw = bsxfun(@times, C_new, 1./aa);
fprintf('Deconvolve and denoise all temporal traces again...\n');
if deconv_flag
    obj.C = deconvTemporal(obj);
else
    obj.C_raw = bsxfun(@minus, obj.C_raw, min(obj.C_raw,[],2)); 
    obj.C = obj.C_raw; 
end
fprintf('Done!\n');

%% upadte b0
if strcmpi(bg_model, 'ring')
    fprintf('Update the constant baselines for all pixels..\n');
    obj.b0_new = obj.P.Ymean-CNMFE_reshape(obj, obj.A*mean(obj.C,2), 2); 
    obj.b0_new = reconstruct_b0(Y_mat(:, :, frame_range(1) : frame_range(2)), obj);
    fprintf('Done!\n');
end

end
function [aa, C_raw] = fast_temporal(Y, A)
%% estimate temporal components using the mean fluorescence
% input:
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components
% output:
%   aa: K*1, energy of each Ai
%   C_raw: K*T, the temporal components before thresholded or being
%   denoised.
% Author: Pengcheng Zhou, Columbia University, 2017
% zhoupc1988@gmail.com

%% options

%% initialization
% roughly initialize C
tmpA = bsxfun(@times, A, 1./max(A, [], 1));
ind_max = bsxfun(@ge, tmpA, 0.5); %max(tmpA, [], 2));
tmp_A = A.*double(ind_max);
aa = sum(tmp_A.^2, 1);
ind = (aa==0);
aa(ind) = inf;
C_raw = bsxfun(@times, tmp_A'*Y, 1./aa');
aa(ind) = 0;
% if any(ind)  % explain the residual using weak neurons
%     tmp_A = A(:, ind);
%     C_raw(ind, :) = (tmp_A'*tmp_A)\(tmp_A'*Y - tmp_A'*A*C_raw);
%     aa(ind) = sum(tmp_A.^2, 1);
% end
end
function b0_ = reconstruct_b0(Y_crop, obj)
            try
                % map data
                mat_data = obj.P.mat_data;
                
                % dimension of data
                d1 = size(Y_crop, 1);
                d2 = size(Y_crop, 2);
                T = size(Y_crop, 3);
                obj.options.d1 = d1;
                obj.options.d2 = d2;
                
                % parameters for patching information
                patch_pos =[[1; d1], [1; d2]];
                % number of patches

            catch
                error('No data file selected');
                b0_= [];
                return;
            end
            
            b0_ = zeros(d1, d2);
            b0_patch = obj.b0;
            tmp_patch = patch_pos;
            r0 = tmp_patch(1);
            r1 = tmp_patch(2);
            c0 = tmp_patch(3);
            c1 = tmp_patch(4);
            b0_(r0:r1, c0:c1) = reshape(b0_patch, r1-r0+1, c1-c0+1);

end













