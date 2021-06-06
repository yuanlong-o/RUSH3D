function [obj, center, Cn, PNR] = initComponents_with_Ain_no_patch(Y_mat, ...
                        obj, Ain, frame_range)
%% initializing spatial/temporal components for calcium imaging data
%% input:
%   Ain:  manually selected ROIs for neurons
%   frame_range: 1 X 2 vector indicating the starting and ending frames
%   save_avi: save the video of initialization procedure
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.
%   use_prev: boolean, use previous initialization or not
%% Output:
%   center: d*2 matrix, centers of all initialized neurons.
%   Cn:     correlation image
%   PNR:    peak to noise ratio

%  modification:
%   only supports
%   no patch involved here

%  Author Yuanlong Zhang. last udpate: 10/15/2020.
%% process parameters
assert(numel(size(Y_mat)) == 3)
try
    % map data
   % dimension of data
    
    d1 = size(Y_mat, 1);
    d2 = size(Y_mat, 2);
    T = size(Y_mat, 3);
    obj.options.d1 = d1;
    obj.options.d2 = d2;
    % frames to be loaded for initialization
    if ~exist('frame_range', 'var')
        frame_range = obj.frame_range;
    end
    if isempty(frame_range)
        frame_range = [1, T];
    else
        frame_range(frame_range<1) = 1;
        frame_range(frame_range>T) = T;
    end
    T = diff(frame_range) + 1;
    obj.frame_range = frame_range;
    outdir = obj.outdir;
    mkdir(obj.outdir)
end

%%
% exporting initialization procedures as a video
% parameter for avoiding using boundaries pixels as seed pixels
options = obj.options;
if ~isfield(options, 'bd') || isempty(options.bd')
    options.bd = options.gSiz;   % boundary pixesl to be ignored during the process of detecting seed pixels
end
bg_ssub = options.bg_ssub;

%% preallocate spaces for saving model variables relating to background components
bg_model = obj.options.background_model;
b = cell(1);
f = cell(1);
% Ymean = cell(nr_patch, nc_patch);
if strcmpi(bg_model, 'ring')
    rr = ceil(obj.options.ring_radius/bg_ssub);    % radius of the ring
    [r_shift, c_shift] = get_nhood(rr, obj.options.num_neighbors);    % shifts used for acquiring the neighboring pixels on the ring

    tmp_patch = [[1; d1], [1; d2]];    % patch position
    tmp_block = [[1; d1], [1; d2]];    % block position
    nr = diff(tmp_patch(1:2)) + 1;
    nc = diff(tmp_patch(3:4)) + 1;
    nr_block = diff(tmp_block(1:2))+1;
    nc_block = diff(tmp_block(3:4))+1;
    b0 = zeros(nr*nc, 1);

    if bg_ssub==1
        [csub, rsub] = meshgrid(tmp_patch(3):tmp_patch(4), tmp_patch(1):tmp_patch(2));
        csub = reshape(csub, [], 1);
        rsub = reshape(rsub, [], 1);
        ii = repmat((1:numel(csub))', [1, length(r_shift)]);
        csub = bsxfun(@plus, csub, c_shift);
        rsub = bsxfun(@plus, rsub, r_shift);
        ind = and(and(csub>=1, csub<=d2), and(rsub>=1, rsub<=d1));
        jj = (csub-tmp_block(3)) * (diff(tmp_block(1:2))+1) + (rsub-tmp_block(1)+1);

        temp = sparse(ii(ind), jj(ind), 1, nr*nc, nr_block*nc_block);
        W = bsxfun(@times, temp, 1./sum(temp, 2));
    else
        d1s = ceil(nr_block/bg_ssub);
        d2s = ceil(nc_block/bg_ssub);

        [csub, rsub] = meshgrid(1:d2s, 1:d1s);
        csub = reshape(csub, [], 1);
        rsub = reshape(rsub, [], 1);
        ii = repmat((1:numel(csub))', [1, length(r_shift)]);
        csub = bsxfun(@plus, csub, c_shift);
        rsub = bsxfun(@plus, rsub, r_shift);
        jj = (csub-1) * d1s + rsub;
        % remove neighbors that are out of boundary
        ind = and(and(csub>=1, csub<=d2s), and(rsub>=1, rsub<=d1s));
        temp = sparse(ii(ind), jj(ind), 1, d1s*d2s, d1s*d2s);
        W = bsxfun(@times, temp, 1./sum(temp, 2)); % this is quite huge!
    end
else
    error('only ring mode supported!')
end

obj.W = W;
obj.b0 = b0;
obj.b = b;
obj.f = f;
clear W b0 b f;    % remove these variables for saving RAM space

%% start initialization
% PNR
PNR = zeros(d1, d2);
% concatenate Y file
% Y_mat = zeros(d1, d2, T);
% for mpatch=1:(nr_patch*nc_patch) %#ok<*UNRCH>
% 	tmp_patch = patch_pos{mpatch};
% 	tmp_block = block_pos{mpatch};
%     Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
% 	Y_mat(tmp_patch(1) : tmp_patch(2), tmp_patch(3) : tmp_patch(4), frame_range(1) : frame_range(2)) =...
%                     Ypatch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1,...
%                     (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1, :);
% end
% Y_mat = double(Y_mat);
% Cn
SI.bg_iter = 10;
SI.pixelwise_window_1 = 10;
SI.pixelwise_window_2 = 200;
SI.pixelwise_poly_index = 2;
SI.Nnum = 15;
SI.detrend_mode = 'pixelwise';
SI.outdir = outdir;
% use imadjust
[Cn, bg_spatial_init,bg_temporal_init] = ...
    summary_image_generation(reshape(Y_mat(:, :, frame_range(1) : frame_range(2)), d1 * d2, T), [d1, d2], SI);

Cn = imadjust(Cn);

% remove backgrounds
Y_mat_wo_bg = reshape(reshape(Y_mat(:, :, frame_range(1) : frame_range(2)), [], T)- bg_spatial_init(:) * bg_temporal_init(:).', d1, d2, T);
% Cin_raw and Cin
Cin_raw = zeros( size(Ain, 2), T);
Cin = zeros( size(Ain, 2), T);
Sin = zeros( size(Ain, 2), T);

for i = 1 : size(Ain, 2)
    curr_A = reshape(full(Ain(:, i)), d1, d2);
    
    % find top and 
    ind_A = find(curr_A);
    [top_left(1), top_left(2)] = ind2sub([d1, d2], ind_A(1));
    [bottom_right(1), bottom_right(2)] = ind2sub([d1, d2], ind_A(end));
    curr_center = [(top_left(1) + bottom_right(1))/2, (top_left(2) + bottom_right(2))/2];
    for j = 1 : T
        curr_Cin_raw(j) = mean(Y_mat_wo_bg(:, :, j) .* curr_A, 'all'); % run a raw background estimation frame by frame
    end
    % deconvolution
    if options.deconv_flag
        % deconv the temporal trace
        [ci, si, deconv_options] = deconvolveCa(curr_Cin_raw, options.deconv_options);  % sn is 1 because i normalized c_raw already
        % save this initialization
        kernel_pars{i} = reshape(deconv_options.pars, 1, []);
        Cin(i, :) = ci;
        Sin(i, :) = si;
        Cin_raw(i, :) = curr_Cin_raw-deconv_options.b;
    else
        Cin(i, :) = curr_Cin_raw;
        Cin_raw(i, :) = curr_Cin_raw;
    end
    
    center(i, :) = curr_center;
end

% Ymean
Ymean = mean(Y_mat(:, :, frame_range(1) : frame_range(2)), 3);

% save the log infomation
log_data.options_0=options;
obj.P.k_options = 1;
obj.P.k_neurons = 0;


%% export the results
obj.A = Ain;
obj.C = Cin;
obj.C_raw = Cin_raw;
if options.deconv_flag
    kernel_pars = cell2mat(reshape(kernel_pars, [], 1));
    obj.S = sparse(Sin);
    obj.P.kernel_pars = kernel_pars;
else
    obj.S = zeros(size(obj.C));
end
obj.Cn = Cn;
K = size(obj.A, 2);
obj.P.k_ids = K;
obj.ids = (1:K);
obj.tags = zeros(K,1, 'like', uint16(0));
obj.P.Ymean = Ymean;
