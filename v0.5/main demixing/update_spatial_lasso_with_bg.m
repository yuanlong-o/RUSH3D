function [A_out, A_bg_out, bg_spatial_out] = update_spatial_lasso_with_bg(Y, ...
    A, A_bg, C, C_bg, bias, bg_spatial, bg_temporal, IND, IND_bg, movie_size, maxIter, sn, q)

%% update spatial components using constrained non-negative lasso with warm started HALS 
%  with consideration about the 

% input: 
%       Y:    d x T,  fluorescence data
%       A:    K element cell,  spatial components
%       C:    K element cell,  temporal components
% bg_spatial: d x 1,  spatial background
% bg_temporal: 1 x T, temporal backgrond
%     IND:    K element cell,  spatial extent for each component
%      sn:    d x 1,  noise std for each pixel
%       q:    scalar, control probability for FDR (default: 0.75), FDR for false discover rate?
% maxIter:    maximum HALS iteration (default: 40)
% options:    options structure

% output: 
%   A_out: d*K, updated spatial components. Note taht K may contain the
%   background term

% last update: 4/23/2020. YZ
%% options for HALS
%norm_C_flag = false;
tol = 1e-3;
repeat = 1;

% penalty related
if nargin < 14 || isempty(q); q = 0.75; end

if nargin<13 || isempty(sn)
    % nosie option
    sn_options.noise_range = [0.25,0.7];
    sn_options.noise_method = 'logmexp';
    sn_options.block_size = [64,64];
    sn_options.split_data = false;
    sn_options.max_timesteps = 3000;
    
    sn = get_noise_fft(Y,sn_options);  
end % sn is the average of power spectrum of Y along time domain

if nargin < 12 || isempty(maxIter); maxIter = 40; end

%% cell to matrix
size_h = movie_size(1);
size_w = movie_size(2);

if iscell(A)
    A_mat = zeros(size(Y, 1), length(A));
    for i = 1 : length(A)
        buf = zeros(size_h, size_w);
        buf_bias = bias{i};
        buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2)) = gather(A{i});
        A_mat(:, i) = buf(:);
    end
end

if iscell(A)
    A_bg_mat = zeros(size(Y, 1), length(A));
    for i = 1 : length(A)
        buf = zeros(size_h, size_w);
        buf_bias = bias{i};
        buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2)) = gather(A_bg{i});
        A_bg_mat(:, i) = buf(:);
    end
end

if iscell(IND)
    IND_mat = false(size(Y, 1), length(IND));
    for i = 1 : length(A)
        buf = false(size_h, size_w);
        buf_bias = bias{i};
        buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2)) = gather(IND{i});
        IND_mat(:, i) = buf(:);
    end
end

if iscell(IND_bg)
    IND_bg_mat = false(size(Y, 1), length(IND));
    for i = 1 : length(A)
        buf = false(size_h, size_w);
        buf_bias = bias{i};
        buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2)) = gather(IND_bg{i});
        IND_bg_mat(:, i) = buf(:);
    end
end

if iscell(C)
    C_mat = zeros(length(C), size(Y, 2));
    for i = 1 : length(C)
        C_mat(i, :) = C{i};
    end
end

if iscell(C_bg)
    C_bg_mat = zeros(length(C_bg), size(Y, 2));
    for i = 1 : length(C_bg)
        C_bg_mat(i, :) = C_bg{i};
    end
end

[d, nr] = size(A_mat);
K = nr + size(bg_spatial, 2) + nr; 
% here K is the neuron number + background component number, nr is the neuron number

%% combine
A_mat = [A_mat, bg_spatial, A_bg_mat];
C_mat = [C_mat; bg_temporal; C_bg_mat];

T = size(C_mat,2);
sn = double(sn);

YC = double(mm_fun(C_mat,Y)); % this is weird, why make C times Y

%% initialization 
V = double(C_mat*C_mat'); 
cc = diag(V);   % array of  squared of l2 norm for all components, this is a easy way to calculate the L2 norm

%% updating (neuron by neuron)
clear C_bg_mat A_bg_mat
miter = 0;
while repeat && miter < maxIter % this is the main iteration. all the step is hals
    A_ = A_mat;
    for k=1:K        % for (each neuron + background component)
        if k <= nr
            lam = sqrt(cc(k)); %max(sqrt(cc(tmp_ind))); % get the l2 norm of V
            tmp_ind = IND_mat(:, k); % note: here IND is a sparse matrix, so only locally update
        elseif k == nr + 1
            lam = 0;
            tmp_ind = true(d, 1);
        elseif k >= nr + 2
            lam = sqrt(cc(k)); %max(sqrt(cc(tmp_ind))); % get the l2 norm of V
            tmp_ind = IND_bg_mat(:, k - (nr + 1)); % note: here IND is a sparse matrix, so only locally update            
        end
        LAM = norminv(q)*sn*lam; % norminv is the inverse of the normcdf, to this is automatically assign the penaly coefficient
        ak = max(0, A_mat(tmp_ind, k)+(full(YC(tmp_ind, k)) - LAM(tmp_ind) ...
            - A_mat(tmp_ind,:)*V(:, k))/cc(k)); % this is a standard HALS step. see the paper for more details
        A_mat(tmp_ind, k) = ak; 
    end
    miter = miter + 1;
    repeat = (sqrt(sum((A_mat(:)-A_(:)).^2)/sum(A_(:).^2)) > tol);    % stop critiera:
end
%% update the variable
% neuron
A_out = [];
A_bg_out = [];

for i = 1 : nr
	buf = A_mat(:, i);
    buf = reshape(buf, size_h, size_w);
    buf_bias = bias{i};
    buf_crop = buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2));
	A_out{i} =buf_crop;
    
	buf = A_mat(:, i + nr + 1);
    buf = reshape(buf, size_h, size_w);
    buf_bias = bias{i};
    buf_crop = buf(buf_bias(1, 1) : buf_bias(2, 1), buf_bias(1, 2) : buf_bias(2, 2));    
    A_bg_out{i} = buf_crop;
end
% background
f = C_mat(nr+1,:); % separate f
Yf = full(YC(:,nr+1)); %Y*f';
b = double(max((double(Yf) - (A_mat(:,1:nr)*double(C_mat(1:nr,:)) + ...
    A_mat(:,nr + 2 : end) * double(C_mat(nr + 2: end,:)))*f')/(f*f'),0)); % separate b
bg_spatial_out = b;
