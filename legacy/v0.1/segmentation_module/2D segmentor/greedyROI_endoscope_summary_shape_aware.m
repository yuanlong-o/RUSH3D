function [results, center, Cn] = greedyROI_endoscope_summary_shape_aware(Y, K, ...
        options,debug_on)
%% a greedy method for detecting ROIs and initializing CNMF. in each iteration,
% it searches the one with large (peak-median)/noise level and large local
% correlation. It's the same with greedyROI_corr.m, but with some features
% specialized for endoscope data
%% Input:
%   Y:  d X T matrx, imaging data
%   K:  scalar, maximum number of neurons to be detected.
%   options: struct data of paramters/options
%       d1:     number of rows
%       d2:     number of columns
%       gSiz:   maximum size of a neuron
%       nb:     number of background
%       min_corr: minimum threshold of correlation for segementing neurons
%   sn:     d X 1 vector, noise level of each pixel
%   debug_on: options for showing procedure of detecting neurons
%   save_avi: save the video of initialization procedure. string: save
%   video; true: just play it; false: interactive mode. (the name of this
%      argument is very misleading after several updates of the code. sorry)

%% Output:
%`      results: struct variable with fields {'Ain', 'Cin', 'Sin', 'kernel_pars'}
%           Ain:  d X K' matrix, estimated spatial component
%           Cin:  K'X T matrix, estimated temporal component
%           Sin:  K' X T matrix, inferred spike counts within each frame
%           kernel_pars: K'X1 cell, parameters for the convolution kernel
%           of each neuron
%       center: K' X 2, coordinate of each neuron's center
%       Cn:  d1*d2, correlation image
%       save_avi:  options for saving avi.

%% Author: Pengcheng Zhou, Carnegie Mellon University. zhoupc1988@gmail.com
% the method is an modification of greedyROI method used in Neuron paper of Eftychios
% Pnevmatikakis et.al. https://github.com/epnev/ca_source_extraction/blob/master/utilities/greedyROI2d.m
% In each iteration of peeling off neurons, it searchs the one with maximum
% value of (max-median)/noise * Cn, which achieves a balance of SNR and
% local correlation.

%  modification:
%      delete all temporal variables
%      set cn (correlation image) to by the input summary image.
%      
%% modified by YZ. last update: 8/28/2020


%% use correlation to initialize NMF
%% parameters
d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSig = options.gSig;    % width of the gaussian kernel approximating one neuron
gSiz = options.gSiz;    % average size of neurons
ssub = options.ssub;
min_threshold = options.min_threshold ;
if (ssub~=1) 
    d1_raw = d1;
    d2_raw = d2;
    T_raw = size(Y, ndims(Y));
    Y = dsData(reshape(double(Y), d1_raw, d2_raw, T_raw), options);
    [d1, d2, ~] = size(Y);
    options.d1 = d1;
    options.d2 = d2;
    gSig = gSig/ssub;
    gSiz = round(gSiz/ssub);
    options.min_pixel = options.min_pixel/(ssub^2);
    options.bd = round(options.bd/ssub);
    
end
false_threshold = options.false_threshold;
min_v_search = options.min_v_search; % serve as a hard threshold
seed_method = options.seed_method; % methods for selecting seed pixels
local_constrast_th = options.local_constrast_th;
% kernel_0 = options.kernel;
min_pixel = options.min_pixel;  % minimum number of pixels to be a neuron

% smin = options.smin;
% boudnary to avoid for detecting seed pixels
try
    bd = options.bd;
catch
    bd = round(gSiz/2);
end

% exporting initialization procedures as a video

% debug mode and exporting results
if ~exist('debug_on', 'var')
    debug_on = false;
end

if ~ismatrix(Y); Y = reshape(Y, d1*d2, []); end;  % convert the 3D movie to a matrix
Y(isnan(Y)) = 0;    % remove nan values
Y = double(Y);
T = size(Y, 2);

%% preprocessing data
% create a spatial filter for removing background
if gSig>0
    if options.center_psf
        psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
        ind_nonzero = (psf(:)>=max(psf(:,1))); % a round shape
        psf = psf-mean(psf(ind_nonzero)); % important here!
%         psf(~ind_nonzero) = 0;
        psf = psf / sum(psf(:));
    else
        psf = fspecial('gaussian', round(gSiz), gSig); %figure, imagesc(psf); axis off, axis equal
    end
    
    % template for background
    local_temp = ones(gSiz, gSiz);
    local_temp(round(gSiz / 2) - floor(gSig) : round(gSiz / 2)+ floor(gSig), ...
               round(gSiz / 2) - floor(gSig) : round(gSiz / 2)+ floor(gSig)) = 0;
    local_temp = local_temp / sum(local_temp(:));
else
    psf = [];
end

% filter the data
if isempty(psf)
    % no filtering
    HY = Y;
else
    HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate'); %figure, imshow(HY, [])
    HY_local = imfilter(reshape(Y, d1,d2,[]), local_temp, 'replicate'); %figure, imshow(HY, [])
end

HY = reshape(HY, d1*d2, []);
% HY_med = median(HY, 2);
% HY_max = max(HY, [], 2)-HY_med;    % maximum projection


% estimate noise level and thrshold diff(HY)

% screen seeding pixels as center of the neuron
Cn = reshape(HY, d1, d2); % correlation image
mask = Cn  > (max(Cn(:)) * 0.05);
mask2 = Cn ./ HY_local > local_constrast_th;
v_search = Cn ./ HY_local .* mask .* mask2; % update base. note v_search shall be 2D
% v_search = Cn .* mask;
ind_search = false(d1*d2,1);  % showing whether this pixel has been searched before
ind_search(v_search==0) = true; % ignore pixels with small correlations or low peak-noise-ratio

figure, imagesc(v_search)
% ignore boundaries pixels when determinging seed pixels
if length(bd) ==1
    bd = ones(1,4)*bd;
end

% cancel the boundaries
ind_bd = false(size(v_search));
ind_bd(1:bd(1), :) = true; % set all boundary to be true. 
ind_bd((end-bd(2)+1):end, :) = true;
ind_bd(:, 1:bd(3)) = true;
ind_bd(:, (end-bd(4)+1):end) = true;

%% start initialization
if ~exist('K', 'var')||isempty(K)
    K = floor(sum(v_search(:)>0)/10);
else
    K = min(floor(sum(v_search(:)>0)/10), K);
end

% initialization
Ain = zeros(d1*d2, K);  % spatial components
center = zeros(K, 2);   % center of the initialized components

if debug_on
    figure('position', [100, 100, 1200, 800], 'color', [1,1,1]*0.9); %#ok<*UNRCH>
    set(gcf, 'defaultAxesFontSize', 20);
    ax_cn = axes('position', [0.04, 0.5, 0.3, 0.4]);
    ax_pnr_cn = axes('position', [0.36, 0.5, 0.3, 0.4]);
    ax_cn_box = axes('position', [0.68, 0.54, 0.24, 0.32]);
    ax_ellipse_fitting = axes('position', [0.68, 0.04, 0.24, 0.32]);
    axes(ax_cn);
    imagesc(Cn);
    %     imagesc(Cn.*PNR, quantile(Cn(:).*PNR(:), [0.5, 0.99]));
    axis equal off; hold on;
    axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
    %     title('Cn * PNR');
    title('Cn');    
end

%% do initialization in a greedy way
searching_flag = true;
k = 0;      %number of found components
[ii, jj] = meshgrid(1:d2, 1:d1);
ind_false_time = 0;
pixel_v = (ii*10+jj)*(1e-10);
while searching_flag && K>0 && ind_false_time < false_threshold
    %% find local maximum as initialization point
    %find all local maximum as initialization point
    
    tmp_d = max(3, round(gSiz/4));
    
    v_search = medfilt2(v_search,3*[1, 1])+pixel_v; % add an extra value to avoid repeated seed pixels within one ROI.
%     figure, imshow(v_search, [])
    v_search(ind_search) = 0; % cancel the last seach
    
    v_max = ordfilt2(v_search, tmp_d^2, true(tmp_d));%figure, imshow(v_max)
    % set boundary to be 0
    v_search(ind_bd) = 0; % cancel the boundary
    
    if strcmpi(seed_method, 'manual') %manually select seed pixels
        tmp_fig = figure('position', [200, 200, 1024, 412]);
        subplot(121); cla;
        imagesc(PNR);  hold on;
        title('PNR');
        plot(center(1:k, 2), center(1:k, 1), '*r');
        axis equal off;
        axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
        colorbar;
        
        subplot(122);
        imagesc(Cn, [min(Cn(:)), 1]); %, [0, max(max(min_v_search(:)*0.99), min_v_search)]);
        colorbar;
        hold on;
        axis equal;
        axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
        drawnow;
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        title('click neuron centers for initialziation');
        xlabel('click invalid pixels to stop', 'color', 'r');
        ind_localmax = zeros(K,1);
        for tmp_k=1:K
            figure(tmp_fig);
            [tmp_x, tmp_y] = ginput(1);
            tmp_x = round(tmp_x); tmp_y = round(tmp_y);
            if isempty(tmp_x)||or(tmp_x<1, tmp_x>d2) || or(tmp_y<1, tmp_y>d1) ||(v_search(tmp_y, tmp_x)==0)
                break;
            end
            plot(tmp_x, tmp_y, '*r', 'linewidth', 2);
            drawnow();
            ind_localmax(tmp_k) = sub2ind([d1,d2], tmp_y, tmp_x);
        end
        close(tmp_fig);
        
        ind_localmax = ind_localmax(1:(tmp_k-1));
        if isempty(ind_localmax)
            break;
        end
    else
        % automatically select seed pixels
        ind_search(v_search<min_v_search) = true;    % avoid generating new seed pixels after initialization
        ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0)); % local max calculation
       if(isempty(ind_localmax)); break; end
   end
    [~, ind_sort] = sort(v_search(ind_localmax), 'descend');
    ind_localmax = ind_localmax(ind_sort); % local max?
    [r_peak, c_peak] = ind2sub([d1,d2],ind_localmax);
    
    %% try initialization over all local maximums
    for mcell = 1:length(ind_localmax)
        % find the starting point
        ind_p = ind_localmax(mcell);
        %         max_v = max_vs(mcell);
        max_v = v_search(ind_p);
        if mcell==1 % first onbe
            img_clim = [0, max_v];
        end
        ind_search(ind_p) = true; % indicating that this pixel has been searched.
        if max_v<min_v_search % all pixels have been tried for initialization
            continue;
        end
        [r, c]  = ind2sub([d1, d2], ind_p);
        
        % select its neighbours for estimation of ai and ci, the box size is
        %[2*gSiz+1, 2*gSiz+1]
        rsub = max(1, -gSiz+r):min(d1, gSiz+r);
        csub = max(1, -gSiz+c):min(d2, gSiz+c);
        [cind, rind] = meshgrid(csub, rsub);
        [nr, nc] = size(cind);
        ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
        HY_box = HY(ind_nhood, :);      % extract temporal component from HY_box
        Y_box = Y(ind_nhood, :);    % extract spatial component from Y_box
        ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
        
        % neighbouring pixels to update after initialization of one
        % neuron. Note this thing is 2 times larger
        rsub = max(1, -2*gSiz+r):min(d1, 2*gSiz+r);
        csub = max(1, -2*gSiz+c):min(d2, 2*gSiz+c);
        [cind, rind] = meshgrid(csub, rsub);
        ind_nhood_HY = sub2ind([d1, d2], rind(:), cind(:));
        [nr2, nc2] = size(cind);
        
        %% show temporal trace in the center
%         if k == 29
%             disp('1')
%         end
        if debug_on
            axes(ax_pnr_cn); cla;
            imagesc(reshape(v_search, d1, d2), img_clim); % [0, max_v]);
            title(sprintf('neuron %d', k+1));
            axis equal off; hold on;
            axis([bd(3), d2-bd(4), bd(1), d1-bd(2)]);
            plot(c_peak(mcell:end), r_peak(mcell:end), '.r');
            plot(c,r, 'or', 'markerfacecolor', 'r', 'markersize', 10);
            
            axes(ax_cn_box);
            imagesc(reshape(Cn(ind_nhood), nr, nc), [0, 1]);
            axis equal off tight;
            title('correlation imagep portion'); 
        end
        
        %% extract ai, ci
        sz = [nr, nc];
        if options.center_psf
            [ai, ind_success] =  extract_a(HY_box, Y_box, ind_ctr, sz, options.spatial_constraints, min_threshold );
        else
            [ai, ind_success] =  extract_a(HY_box, Y_box, ind_ctr, sz, options.spatial_constraints, min_threshold );
        end
        if any(isnan(ai)); ind_success=false; end
%         if sum(ai)<=min_pixel; ind_success = false; end % interesting, it is pixel intensity level
        %         if max(ci_raw)<min_pnr;
        %             ind_success=false;
        %         end
        if sum(ai(:)>0)<min_pixel; ind_success=false; end
        
        % shape judgement
        ai_2d = reshape(ai, sz(1), sz(2));
        % ellipse fiting
        BW = ai_2d > (min_threshold * (max(ai_2d(:)) - min(ai_2d(:)))+ min(ai_2d(:)));
        B = bwboundaries(BW);
        B_1 = B{1}(:, 1);
        B_2 = B{1}(:, 2);
        
        

%         B = BW - imerode(BW, true(3));
%         B_k = find(B);
%         [B_1, B_2] = ind2sub([sz(1), sz(2)], B_k);
        try       
            if debug_on
                axes(ax_ellipse_fitting);
                imagesc(ai_2d > min_threshold * max(ai_2d(:)), [0, 1]);
                axis equal off tight;
                title('correlation imagep portion'); 
                [ellipse_t,rotated_ellipse] = fit_ellipse( B_1, B_2, ax_ellipse_fitting);
            else
                [ellipse_t] = fit_ellipse( B_1, B_2);
                R = [ cos(ellipse_t.phi), sin(ellipse_t.phi); -sin(ellipse_t.phi), cos(ellipse_t.phi) ];
                theta_r         = linspace(0,2*pi, 500);
                ellipse_x_r     = ellipse_t.X0 + ellipse_t.a*cos( theta_r );
                ellipse_y_r     = ellipse_t.Y0 + ellipse_t.b*sin( theta_r );
                rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
            end
            
            
            BW_ellipse = zeros(sz(1), sz(2));
            BW_ellipse(sub2ind(sz, round(max(min(rotated_ellipse(1, :), sz(1)), 1)), ...
                                   round(max(min(rotated_ellipse(2, :), sz(1)), 1)))) = 1;
            BW_ellipse = imfill(BW_ellipse, 'holes');  % 
            ellipse_area = ellipse_t.a * ellipse_t.b * pi;
            if isempty(ellipse_t.a) || isempty(ellipse_t.b) % fit success
                ind_success= false;
            elseif sum(BW .* BW_ellipse, 'all') < 0.5 * sum(BW, 'all')% fit area overlap with the neuron area
                ind_success= false;
            elseif ellipse_t.a / ellipse_t.b > 5 || ellipse_t.b / ellipse_t.a > 5 % can not have too weird shape
                ind_success= false;
            elseif ellipse_area > 2 * sum((ai>0) .* BW(:), 'all') ||sum((ai>0) .* BW(:), 'all') > ellipse_area * 2% the shape can not be too weird
                ind_success = false;                
            end
        catch
            ind_success = false;
        end

        
        if ind_success
            k = k+1;
            
            Ain(ind_nhood, k) = ai;
            center(k, :) = [r, c];
            
            % avoid searching nearby pixels
            ind_search(ind_nhood(ai>max(ai)*0.5)) = true;
            
            % update the raw data
            Y(ind_nhood) = Y_box - ai; % directly minus
            % update filtered data
            if isempty(psf)
                Hai = reshape(Ain(ind_nhood_HY, k), nr2, nc2);
            else
                Hai = imfilter(reshape(Ain(ind_nhood_HY, k), nr2, nc2), psf, 'replicate');
            end
            HY_box = HY(ind_nhood_HY) - Hai(:);
            %             HY_box = bsxfun(@minus, HY_box, median(HY_box, 2));
            HY(ind_nhood_HY) = HY_box;
            
            % update the maximum projection of HY
            
            % update search value
            v_search(ind_bd) = 0;
            v_search(ind_search) = 0;
            if debug_on
                axes(ax_cn);
                plot(c, r, '.r'); % plot the sensored points
                axes(ax_pnr_cn); % 
                plot(c,r, 'or');
                axes(ax_cn_box);
                imagesc(reshape(ai, nr, nc));
                axis equal off tight;
                title('spatial component');
                drawnow()
            end
        else
            ind_false_time = ind_false_time +1;
            continue;
        end
        
        %% display results

               
        if mod(k, 10)==0
            fprintf('%d neurons have been detected\n', k);
        end
        
        if k==K
            searching_flag = false;
            break;
        end
    end
end
center = center(1:k, :);
results.Ain = sparse(Ain(:, 1:k));


% Cin(Cin<0) = 0;

if ssub~=1
    Ain = reshape(full(results.Ain), d1,d2, []);
    if ~isempty(Ain)
        results.Ain = sparse(reshape(imresize(Ain, [d1_raw, d2_raw]), d1_raw*d2_raw, []));
    else
        results.Ain = sparse(zeros(d1_raw*d2_raw, 0));
    end
    Cn =imresize(Cn, [d1_raw, d2_raw]);
    center = center*ssub-1; 
end


end
