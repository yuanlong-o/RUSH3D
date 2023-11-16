function seed_param = zrange_def(std_volume,seed_param)

z_tmp = floor(size(std_volume,3)/10);
% gradient
grad_std_volume = zeros(size(std_volume));
for z = 1:size(std_volume,3)
    [grad_x,grad_y] = gradient(std_volume(:,:,z));
    grad_std_volume(:,:,z) = grad_x.^2 + grad_y.^2;
end
gradz_energy = squeeze(mean(grad_std_volume,[1,2]));
[gpek,gpek_loc] = findpeaks(gradz_energy(z_tmp : end - z_tmp));
[~,zmax_loc] = max(gradz_energy(z_tmp : end - z_tmp));
if ~isempty(gpek_loc)
    [~,tmp_gmax_loc] = max(gpek(:));
    zgrad_loc = gpek_loc(tmp_gmax_loc);
else
    zgrad_loc = zmax_loc;
end
zgrad_loc = zgrad_loc + z_tmp - 1;
% energy
stdz_energy = squeeze(mean(std_volume,[1,2]));

[zpek,zpek_loc] = findpeaks(stdz_energy(z_tmp : end - z_tmp));
[~,zmax_loc] = max(stdz_energy(z_tmp : end - z_tmp));
if ~isempty(zpek_loc)
    [~,tmp_zmax_loc] = max(zpek(:));
    zene_loc = zpek_loc(tmp_zmax_loc);
else
    zene_loc = zmax_loc;
end
zene_loc = zene_loc + z_tmp - 1;

zmid_loc = round(zgrad_loc * 0.6 + zene_loc * 0.4);
start_ind = max(zmid_loc - floor(100/(seed_param.per_slice_depth*1e6)),3); % vessel - 60um
end_ind = min(zmid_loc + floor(200/(seed_param.per_slice_depth*1e6)),size(std_volume,3)-3); % vessel + 200um
seed_param.start_ind = start_ind;
seed_param.end_ind = end_ind;

end