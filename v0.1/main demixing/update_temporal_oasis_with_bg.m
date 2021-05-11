
function [C_raw_out, C_bg_out, bg_temporal_out] = ...
    update_temporal_oasis_with_bg(Y,A, A_bg, C, C_bg, bias, bg_spatial, bg_temporal, ...
    IND, IND_bg, movie_size, maxIter, sn, q)

%% update the temporal components using constrained HALS and oasis
%  with background bg. deconvolution only work for the bg removed signals

% input: 
%       Y:    d x T,  fluorescence data
%       A:    K element cell,  spatial components
%       A_bg:  rings for each neuron
%       C:    K element cell,  temporal components
% bg_spatial: d x 1,  spatial background
% bg_temporal: 1 x T, temporal backgrond
%     IND:    K element cell,  spatial extent for each component
%bg_mask_init : background mask, e.g. a ring for each neuron. 
%      sn:    d x 1,  noise std for each pixel
%       q:    scalar, control probability for FDR (default: 0.75), FDR for false discover rate?
% maxIter:    maximum HALS iteration (default: 40)
% options:    options structure

% output: 
%   C_raw_out: K*T, updated spatial components. Note taht K may contain the
%   C_bg_out
%   C_denoised_out
%   Spike_mat
%   bg_temporal_out

%   background term

% update: add spike matrix as output

% last update: 4/15/2020. YZ
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
        buf = C{i};
        C_mat(i, :) = buf(:).';
    end
end

if iscell(C_bg)
    C_bg_mat = zeros(length(C_bg), size(Y, 2));
    for i = 1 : length(C_bg)
        buf = C_bg{i};
        C_bg_mat(i, :) = buf(:).';
    end
end

[d, nr] = size(A_mat);
% here K is the neuron number + background component number, nr is the neuron number

%% combine
A_mat = [A_mat, bg_spatial];
C_mat = [C_mat; bg_temporal];

T = size(C_mat,2);
sn = double(sn);

for i = 1 : nr
    buf = A_mat(:, i);
    tmp_ind = IND_mat(:, i);
    buf(~tmp_ind) = 0;
    A_mat(:, i) = buf;
end

for i = 1 : nr
    buf = A_bg_mat(:, i);
    tmp_ind = IND_bg_mat(:, i);
    buf(~tmp_ind) = 0;
    A_bg_mat(:, i) = buf;

end
% figure, imshow(reshape(A_bg_mat(:, 1), [855, 840]), [])
%% HALS
%  HALS for neuron
[C_mat_combine, ~] = HALS_temporal(Y, [A_mat, A_bg_mat], [C_mat; C_bg_mat], maxIter, tol, true, false);


C_mat = C_mat_combine(1 : nr, :);
f = C_mat_combine(nr + 1, :);
C_bg_mat = C_mat_combine(nr + 2 : end, :);

%% OASIS
% C_pure = C_mat - C_bg_mat; % question: what to do if there are negative signals?
% decimate = 1;
% % optimize_b = false;
% % optimize_g = true;
% % maxIter = 100; % which is only been used to get gama and background
% gama = 0.9;
% % lam = 0.05;
% maxIter_oasis = 100;
% Spike_mat = zeros(nr, size(C_mat, 2));
% 
% C_denoised_mat = C_mat;
% 
% 
% % do normalization before OASIS
% for i = 1 : nr
% 
%     y = C_pure(i, :);
%     [c, s, b, g, active_set] = foopsi_oasisAR1(y, gama, lam, false,...
%         true, decimate, maxIter_oasis);
%     
%     % 
% %     figure(102), plot(c)
% %     figure(103), plot(y)
%     C_denoised_mat(i, :) = c;
%     Spike_mat(i, :) = s;
% end
% figure, plot(C_mat(end, :))
% hold on, plot(C_bg_mat(end, :))
% hold on, plot(C_pure(end, :))
%% update the variable
% neuron
C_raw_out = [];
C_bg_out = [];
for i = 1 : nr
    C_raw_out{i} = C_mat(i, :);
    C_bg_out{i} = C_bg_mat(i, :);
end
% background
bg_temporal_out = f;
