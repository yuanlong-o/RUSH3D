function [bg_spatial_init, bg_temporal_init] = initialize_bg(processed_video, bg_iter)
valid_wigner = size(processed_video, 4);
for j = 1 : valid_wigner
    sensor_movie = processed_video(:, :, :, j);
    [~, ~, T] = size(sensor_movie);
    sensor_movie = reshape(sensor_movie, [],T);
    [spatial_buf, temporal_buf] = rank_1_factorization(sensor_movie, bg_iter);
    
    bg_spatial_init(:, j) = spatial_buf;
    bg_temporal_init(:, j) = temporal_buf;
end
bg_temporal_init = bg_temporal_init.';
bg_temporal_init = mean(bg_temporal_init, 1);
end