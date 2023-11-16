function [summary_image, bg_spatial_init,bg_temporal_init] = ...
    summary_image_generation(sensor_movie, movie_size, SI)

%% module for loading sensory image of LFM.
%   Input
%   sensor_movie: xyt data
%   patch_info: size and location info of current patch
%   SI: configuration structure
%       SI.detrend_mode: summary image calculation mode, you can choose
%       local, global, and pixelwise
%       detailed parameter see the parser section
%   Output
%   summary_image: summarized image (xy data)
%   bg_spatial_init,bg_temporal_init: initilaized rank-1 spatial and
%   temporal background activities

%   last update: 5/31/2020. YZ

%% parser
Nnum = SI.Nnum;
size_h = movie_size(1);
size_w = movie_size(2);
% movie_size = [size_h, size_w];
outdir = SI.outdir;
buf = size(sensor_movie);
frameN = buf(end);
% parser with different method
detrend_mode = SI.detrend_mode;
if strcmp(detrend_mode, 'local')
    block_size = SI.block_size;
end

if strcmp(detrend_mode, 'global')
    slide_windows_size = SI.global_detrend_delta;
    frame_step = SI.frames_step;
    bg_iter = SI.bg_iter;
end

if strcmp(detrend_mode, 'pixelwise')  
    bg_iter = SI.bg_iter;
    window_1 = SI.pixelwise_window_1;
    window_2 = SI.pixelwise_window_2;
    poly_index = SI.pixelwise_poly_index;
end
    
%% main body

if strcmp(detrend_mode, 'local')
    [adj_local_std, local_adjust, local_bg_s, local_bg_t, local_bg_block]...
        = local_bg_remove(Nnum, reshape(sensor_movie, size_h, size_w, []), block_size);
    figure, imagesc(cell2mat(adj_local_std)), axis equal, axis off, title('std without local bg')
    summary_image = adj_local_std;
elseif strcmp(detrend_mode,'global')
    % get baseline
    baseline_raw = squeeze(mean(sensor_movie,1))';
    smooth_window_span = 2 *  slide_windows_size / frame_step;
    % smooth baseline
    baseline = smooth(baseline_raw, smooth_window_span, 'sgolay', 3);
    figure; hold on; plot(baseline_raw); plot(baseline); title('Frame means (post bg subtract), raw + trend fit'); hold off;
    print(fullfile(outdir, ['_trend_fit.pdf']), '-dpdf', '-r300');
    
    % devide the base line
    sensor_movie_max = max(sensor_movie(:));
    sensor_movie = sensor_movie/sensor_movie_max;
    
    [bg_spatial_init, bg_temporal_init] = rank_1_factorization(sensor_movie,bg_iter);
    summary_image = compute_std_image(sensor_movie, bg_spatial_init, bg_temporal_init); % also input the spatial and temporal background
    
    bg_spatial_init = reshape(bg_spatial_init, movie_size(1:2)); 
    summary_image = reshape(summary_image, movie_size(1:2));
elseif strcmp(detrend_mode, 'pixelwise')    
    % rank 1 background removal
	sensor_movie_max = max(sensor_movie(:));
    sensor_movie = sensor_movie/sensor_movie_max;
    
    img_stack = reshape(sensor_movie, movie_size(1), movie_size(2), []);
    
    [bg_spatial_init,bg_temporal_init]=rank_1_factorization(sensor_movie, bg_iter);
    bg_spatial_init = reshape(bg_spatial_init, movie_size(1:2)); 
     
    % detrend
%     window_1 = 200; % large window, which is much larger than a calcium transients
    smoothed_value = smoothdata(img_stack, 3,'movmean', window_1);

    detrend_data = img_stack - smoothed_value;
    figure, plot(squeeze(detrend_data(round(end / 2), round(end / 2), :)))

    % do a smooth to avoid noise contribute to std
%     window_2 = 20;    
    smoothed_value2 = smoothdata(detrend_data, 3,'movmean', window_2);
    hold on, plot(squeeze(smoothed_value2(round(end / 2), round(end / 2), :)))
	print(fullfile(SI.outdir, [ '_example_pixelwise_fit.pdf']), '-dpdf', '-r300');
   
    
    tmp = smoothed_value2.^poly_index;
    % std
    summary_image =compute_std_image(reshape(tmp, [], frameN), ...
        zeros(movie_size(1) * movie_size(2), 1), zeros(1, frameN));
    summary_image = reshape(summary_image,  movie_size(1), movie_size(2));   
    clear img_stack detrend_data smoothed_value smoothed_value2 tmp 
end 


%% save and return
summary_image = summary_image / max(summary_image(:));
saveastiff(im2uint16(summary_image), ...
    fullfile(SI.outdir, ['_summary_img.tif']))
figure, imagesc(summary_image), axis equal, axis off, title('summary image without global bg')
saveas(gca, fullfile(outdir, ['_summary_image.png']))

figure, imagesc(bg_spatial_init), axis equal, axis off, title('rank-1 spatial bg')
saveas(gca, fullfile(outdir, ['_bg_spatial_rank1.png']))

figure, plot(bg_temporal_init), title('rank-1 temporal bg')
saveas(gca, fullfile(outdir, ['_bg_temporal_rank1.png']))
