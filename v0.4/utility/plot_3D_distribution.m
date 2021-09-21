function plot_3D_distribution(valid_seg, movie_size, x1objspace, x3objspace, sz)
%% this funtion is used to save the 3D neuron location in the volume.
%  last upate: 4/1/2020. YZ

neuron_position = zeros(length(valid_seg), 3);
x_pixel_size = abs(x1objspace(2) - x1objspace(1)) * 1e6;
z_pixel_size = abs(x3objspace(2) - x3objspace(1))* 1e6;

FOV_x = x_pixel_size * movie_size(1);
FOV_y = x_pixel_size * movie_size(2);
FOV_z = (max(x3objspace) - min(x3objspace)) * 1e6;

for i = 1 : length(valid_seg)
   pos = valid_seg{i, 2};
   pos_mean = mean(pos, 1);
   pos_mean(1) = pos_mean(1) * x_pixel_size;
   pos_mean(2) = pos_mean(2) * x_pixel_size;
   pos_mean(3) = - x3objspace(round(pos_mean(3))) * 1e6;
   neuron_position(i, :) = pos_mean;
end
scatter3(neuron_position(:, 1) , ...
        neuron_position(:, 2),...
        neuron_position(:, 3), sz, 'filled')
ax = gca;
ax.BoxStyle = 'full';
box on
axis equal
axis vis3d % important, not change the size
xlabel('x axis [um]')
ylabel('y axis [um]')
zlabel('z axis [um]')

% xticks([1, ceil(FOV_x /2), FOV_x])
% xticklabels({'0', ...
%     sprintf('%.0f', FOV_x / 2),...
%     sprintf('%.0f', FOV_x)})
% yticks([1, ceil(FOV_y /2), FOV_y])
% yticklabels({'0', ...
%     sprintf('%.0f', FOV_y / 2 ),...
%     sprintf('%.0f', FOV_y)})
% zticks([1, ceil(FOV_z /2), FOV_z])
% zticklabels({'0', ...
%     sprintf('%.0f', FOV_z / 2 ),...
%     sprintf('%.0f', FOV_z)})

end