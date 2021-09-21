function [r1_shift, r2_shift, r3_shift] = get_nhood_3D(radius, k, thickness)
% find the pixel locations
%% determine neibours of each pixel
if ~exist('k', 'var')
    k = [];
end
r1sub = (-radius):(radius);      % row subscript
r2sub = r1sub;      % column subscript
r3sub = r1sub;

[r2ind, r1ind, r3ind] = meshgrid(r2sub, r1sub, r3sub);
R = sqrt(r2ind.^2+r1ind.^2 +r3ind.^2);
neigh_kernel = (R>=radius) .* (R<radius+thickness);  % kernel representing the selected neighbors

[r1_shift, r2_shift, r3_shift] = ind2sub(size(neigh_kernel),find(neigh_kernel == 1));
r1_shift = reshape(r1_shift - radius - 1, 1, []);
r2_shift = reshape(r2_shift - radius - 1, 1, []); % since in the center
r3_shift = reshape(r3_shift - radius - 1, 1, []);

if isempty(k) || (k>length(r1_shift))
    return;
else
    temp = atan2(r1_shift, r2_shift);
    [~, ids] = sort(temp);
    ind = round(linspace(1, length(ids), k));
    r1_shift = r1_shift(ids(ind));
    r2_shift = r2_shift(ids(ind));
    r3_shift = r3_shift(ids(ind));
end