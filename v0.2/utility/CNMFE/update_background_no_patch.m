function obj = update_background_no_patch(Y_mat, obj)
%% update the background related variables in CNMF framework
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%  cancel mpatch and all patch loop
%% modified by YZ. last update: 10/15/2020.
%% process parameters

try   
    % dimension of data
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
fprintf('\n-----------------UPDATE BACKGROUND---------------------------\n');

% frames to be loaded for initialization
frame_range = obj.frame_range;
T = diff(frame_range) + 1;

% threshold for detecting large residuals
thresh_outlier = obj.options.thresh_outlier;

% options
options = obj.options;
nb = options.nb;
bg_ssub = options.bg_ssub;
bg_model = options.background_model;
with_projection = options.bg_acceleration;

% previous estimation
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;


%% check whether the bg_ssub was changed
if strcmpi(bg_model, 'ring')
    tmp_block = block_pos;    % block position
    nr_block = diff(tmp_block(1:2))+1;
    nc_block = diff(tmp_block(3:4))+1;
    [~, temp] = size(W);
    [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
    if temp~=d1s*d2s
        rr = ceil(obj.options.ring_radius/bg_ssub);    % radius of the ring
        [r_shift, c_shift] = get_nhood(rr, obj.options.num_neighbors);    % shifts used for acquiring the neighboring pixels on the ring

        tmp_patch = patch_pos;    % patch position
        tmp_block = block_pos;    % block position
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
            W = bsxfun(@times, temp, 1./sum(temp, 2));
        end

    end
end

%% start updating the background

tmp_block = block_pos;

% find the neurons that are within the block
mask = zeros(d1, d2);
mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;

ind = (reshape(mask(:), 1, [])* obj.A>0);
A= obj.A(logical(mask), ind);
C = obj.C(ind, :);
temp = obj.P.sn(logical(mask));
if bg_ssub==1
    sn = temp;
else
    nr_block = diff(tmp_block(1:2))+1;
    nc_block = diff(tmp_block(3:4))+1;
    sn= imresize(reshape(temp, nr_block, nc_block), 1/bg_ssub, 'nearest')*bg_ssub;
end


% check whether this is the first run of updating background components
if strcmpi(bg_model, 'ring')
    flag_first = (length(unique(W(1, :)))==2);
else
    flag_first = (mean2(b)==0);
end



tmp_patch = patch_pos;
tmp_block = block_pos;

A_block = A;
sn_block = sn; % sn for spike
C_block = C;

% stop the updating B because A&C doesn't change in this area

% use ind_patch to indicate pixels within the patch and only
% update (W, b0) corresponding to these pixels
ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;

% pull data
Ypatch = Y_mat(:, :, frame_range(1) : frame_range(2));
[nr_block, nc_block, T_block] = size(Ypatch);
if strcmpi(bg_model, 'ring')
    % get the previous estimation
    W_old = W;
    Ypatch = reshape(Ypatch, [], T_block);

    % run regression to get A, C, and W, b0
    if bg_ssub==1
        sn_patch = sn_block(ind_patch);
        [W, b0] = fit_ring_model(Ypatch, A_block, C_block, W_old, thresh_outlier, sn_patch, ind_patch, with_projection);
    else
        % downsapmle data first
        temp = reshape(double(Ypatch)-A_block*C_block, nr_block, nc_block, T_block);
        tmp_b0 = mean(temp, 3);
        b0 = tmp_b0(ind_patch);
        Ypatch = imresize(temp, 1./bg_ssub, 'nearest');
        Ypatch = reshape(Ypatch, [], T_block);

        [W, ~] = fit_ring_model(Ypatch, [], [], W_old, thresh_outlier, sn_block(:), [],  with_projection);
        %                 tmp_b0 = imresize(reshape(tmp_b0, size(sn_block)), [nr_block, nc_block]);
        %                 b0{mpatch} = tmp_b0(ind_patch(:));
    end
else
    error('only ring mode supported!')
end

obj.b = b;
obj.f = f;
obj.b0 = b0;
obj.W = W;
obj.b0_new = reconstruct_b0(Y_mat(:, :, frame_range(1) : frame_range(2)), obj);
obj.A_prev = obj.A;
obj.C_prev = obj.C;

%% save the results to log
fprintf('Finished updating background using %s model.\n', bg_model);
end

function b0_ = reconstruct_b0(Y_crop, obj)
            try
                
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
