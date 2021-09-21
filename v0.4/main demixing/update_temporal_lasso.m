function [C_mat, C_bg_mat , bg_temporal_out] = ...
    update_temporal_lasso(processed_video,A_mat, A_bg_mat, C, C_bg, bg_spatial, bg_temporal, ...
    IND_mat, IND_bg_mat, maxIter, sn, q)

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

movie_size = size(processed_video);
%% concatenate different views
size_h = movie_size(1);
size_w = movie_size(2);
Y = zeros(size_h * size_w * movie_size(4), movie_size(3), 'single');
for i = 1 : movie_size(4)
    buf = processed_video(:, :, :, i);
    buf = reshape(buf, [], movie_size(3));
    Y(size_h * size_w  * (i - 1) + 1 : size_h * size_w  * i, :) = buf;
end
%% options for HALS
%norm_C_flag = false;
tol = 1e-3;
repeat = 1;

% penalty related
if nargin < 12 || isempty(q); q = 0.75; end

if nargin<11 || isempty(sn)
    % nosie option
    sn_options.noise_range = [0.25,0.7];
    sn_options.noise_method = 'logmexp';
    sn_options.block_size = [64,64];
    sn_options.split_data = false;
    sn_options.max_timesteps = 3000;
    
    sn = get_noise_fft(Y,sn_options);  
end % sn is the average of power spectrum of Y along time domain

if nargin < 10 || isempty(maxIter); maxIter = 40; end

%% cell to matrix

% process A
% A_mat = cell_process_A(A, movie_size);
% clear A
% A_bg_mat = cell_process_A(A_bg, movie_size);
% clear A_bg
% IND_mat = cell_process_A(IND, movie_size);
% IND_bg_mat = cell_process_A(IND_bg, movie_size);

% process bg
bg_spatial = reshape(bg_spatial, [], 1);

[d, nr] = size(A_mat);
K = nr + size(bg_spatial, 2) + nr; 
% here K is the neuron number + background component number, nr is the neuron number

%% combine
A_mat = [A_mat, bg_spatial];
C_mat = [C; bg_temporal];

T = size(C_mat,2);
sn = double(sn);

% figure, imshow(reshape(A_bg_mat(:, 1), [855, 840]), [])
%% HALS
%  HALS for neuron
[C_mat_combine, ~] = HALS_temporal_ext(Y, full([A_mat, A_bg_mat]), full([C_mat; C_bg]), maxIter, tol, true, false);


C_mat = C_mat_combine(1 : nr, :);
f = C_mat_combine(nr + 1, :);
C_bg_mat = C_mat_combine(nr + 2 : end, :);


%% update the variable
% neuron
% background
bg_temporal_out = f;
