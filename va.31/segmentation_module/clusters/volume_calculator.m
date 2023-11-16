function [dia_volume, volume_size] = volume_calculator(centers, patch_size, size_h, size_w, size_z)
    % calculate volume based on the dilate
    volume = false(size_h, size_w, size_z);
    ind = sub2ind([size_h, size_w, size_z], centers(:, 1), centers(:, 2), centers(:, 3));
    volume(ind) = 1;
    % build template
    se = strel('cuboid', [patch_size, patch_size, 1]);
    dia_volume = imdilate(volume, se); % segments
    volume_size = sum(dia_volume, 'all');
end