function [S_init, S_shell_init,S_mask_init,  S_shell_mask_init] = seed_generation_module(psf, valid_seg, ...
                                                                volume_size, wdf_size, specify_wigner, shell_radius, seed_param)
%% this file will generate an iteration seed for each neuron.
%  input 3D valid seg valid_seg
%  

%  last update: 5/23/2021. YZ


%% parameters
vol_ds = seed_param.vol_ds;
wdf_ds = seed_param.wdf_ds;
NN_ds = seed_param.NN_ds;

%% initialization

N_seg = size(valid_seg, 1);

S_init = cell(N_seg, size(specify_wigner, 1));
S_mask_init = cell(N_seg, size(specify_wigner, 1));

S_shell_init = cell(N_seg, size(specify_wigner, 1));
S_shell_mask_init = cell(N_seg, size(specify_wigner, 1));

ds_ratio = NN_ds * wdf_ds / vol_ds;

% psf_d = imresize(psf, 1/ds_ratio);
psf_d = imresize(psf, [floor(size(psf,1)/ds_ratio/2)*2+1,floor(size(psf,2)/ds_ratio/2)*2+1],'cubic');
discard_array = zeros(size(valid_seg, 1), 1);
%% forward projection
for j = 1 : size(valid_seg, 1)
    if mod(j, 10) == 0
        fprintf('%d of %d processed \n', j, size(valid_seg, 1))
        
    end
    % register the valid_seg in 3D space
    curr_seg = valid_seg{j, 1}; % still a cell array
    curr_pos = valid_seg{j, 2}; % a matrix
    num_comp_in_seg = length(curr_seg);
    patch_size = size(curr_seg{1}, 1);   
    
    
    % restore the seg in 3D
    [A_in, center_A_in] = generate_Ain_3D(curr_seg, curr_pos , volume_size(1), volume_size(2), volume_size(3));
    
    % generate a shell
    A_in_shell = build_shell_3D(center_A_in , volume_size, shell_radius, ds_ratio);
    
    % down sampling
    A_in = imresize(A_in,[wdf_size(1),wdf_size(2)],'nearest');
    %A_in = imresize(A_in, 1 / ds_ratio, 'nearest');
    A_in_shell = imresize(A_in_shell,[wdf_size(1),wdf_size(2)],'nearest');
    %A_in_shell = imresize(A_in_shell, 1 / ds_ratio, 'nearest');
        
    for i = 1 : size(specify_wigner, 1)
        curr_vind = specify_wigner(i, 1);
        % sample
        HXguess = prop_to_target_wigner(A_in, psf_d, curr_vind);
        HXguess_shell = prop_to_target_wigner(A_in_shell, psf_d, curr_vind);
        
        HXguess_mask = prop_to_target_wigner(single(A_in > 0.1 * max(A_in(:))), psf_d, curr_vind);
        HXguess_mask = HXguess_mask > 0.1 * max(HXguess_mask(:));
        HXguess_mask_shell = HXguess_shell;
        HXguess_mask_shell = HXguess_mask_shell > 0.05 * max(HXguess_mask_shell(:));
        
        % record
        S_init{j, i} = sparse(double(HXguess)); % 2d here
        S_shell_init{j, i} = sparse(double(HXguess_shell));
        
        S_mask_init{j, i} = sparse(double(HXguess_mask));
        S_shell_mask_init{j, i} = sparse(double(HXguess_mask_shell));
        
        % safety check
        if max(HXguess(:)) == 0
            discard_array(j) = 1;
        end
        
    end
end
%% discard
discard_array = find( discard_array);
S_init(discard_array, :) = [];
S_shell_init(discard_array, :) = []; 
S_mask_init(discard_array, :) = [];
S_shell_mask_init(discard_array, :) =[];

end

%% utility function
function HXguess = prop_to_target_wigner(Xguess, psf_z, curr_v)

    [size_h, size_w, ~] = size(Xguess);
    tmp = gpuArray(psf_z(:,:,curr_v,:));%%uv?psf???phase space
    tmp = squeeze(tmp);
    [size_psf_h, size_psf_w, ~] = size(tmp);
    size_h_L = max(size_h, size_psf_h);
    size_w_L = max(size_w, size_psf_w);
    % pad to large psf
    large_psf = zeros(size_h_L, size_w_L, size(tmp, 3));
    large_Xguess = zeros(size_h_L, size_w_L, size(tmp, 3));
    large_psf(floor((size_h_L + 1) / 2) - (size_psf_h - 1) / 2 : floor((size_h_L + 1) / 2) + (size_psf_h  - 1) / 2, ...
        floor((size_w_L + 1) / 2) - (size_psf_w - 1) / 2 : floor((size_w_L + 1) / 2) + (size_psf_w - 1) / 2, :) = tmp;
    large_Xguess(floor((size_h_L + 1) / 2) - floor((size_h - 1) / 2) : floor((size_h_L + 1) / 2) - floor((size_h - 1) / 2) + size_h -1, ...
        floor((size_w_L + 1) / 2) - floor((size_w - 1) / 2) : floor((size_w_L + 1) / 2) - floor((size_w - 1) / 2) + size_w -1, :) = Xguess;
    FFTPSF = fftn(flip(large_psf,3));
    tmptmp = fftshift(fftshift(large_Xguess,1),2);%%Xguess????volume
    HXguess =sum((fftn(tmptmp).*FFTPSF),3)./size(psf_z,4);
    HXguess=abs(ifftn(HXguess)); % in large size.
    HXguess = HXguess(floor((size_h_L + 1) / 2) - floor((size_h - 1) / 2) : floor((size_h_L + 1) / 2) - floor((size_h - 1) / 2) + size_h -1, ...
        floor((size_w_L + 1) / 2) - floor((size_w - 1) / 2) : floor((size_w_L + 1) / 2) - floor((size_w - 1) / 2) + size_w -1);
end

function [A_in, center_A_in] = generate_Ain_3D(seg, center, size_h, size_w, size_z)
    A_in = zeros(size_h, size_w, size_z);
    for i = 1 : length(seg)
        curr_center = center(i, 1 : 3);
        [curr_patch_size_h, curr_patch_size_w] = size(seg{i});
        curr_patch_top_left = [curr_center(1) - floor(curr_patch_size_h / 2), ...
                               curr_center(2) - floor(curr_patch_size_w / 2)];
        curr_patch_top_left(curr_patch_top_left < 1) = 1;

        curr_patch_bottom_right = [curr_patch_top_left(1) + curr_patch_size_h - 1, ...
                               curr_patch_top_left(2) + curr_patch_size_w - 1];
        
        if curr_patch_bottom_right(1) > size_h;  curr_patch_bottom_right(1) = size_h; end
        if curr_patch_bottom_right(2) > size_w;  curr_patch_bottom_right(2) = size_w; end
            
        update_patch_size_h = curr_patch_bottom_right(1) - curr_patch_top_left(1) + 1;
        update_patch_size_w = curr_patch_bottom_right(2) - curr_patch_top_left(2) + 1;
        
        A_in(curr_patch_top_left(1) : curr_patch_bottom_right(1),...
             curr_patch_top_left(2) : curr_patch_bottom_right(2),...
             curr_center(3)) = ...
        A_in(curr_patch_top_left(1) : curr_patch_bottom_right(1),...
             curr_patch_top_left(2) : curr_patch_bottom_right(2),...
             curr_center(3)) + ...
             seg{i}(1 : update_patch_size_h, 1 : update_patch_size_w);
         
         
    end
    % average the patch
    A_in = A_in / length(seg);
    center_A_in = round(mean( center, 1));
end


function buf2 = build_shell_3D(center_A_in , movie_size, shell_radius, thickness)
	[r1_shift, r2_shift, r3_shift] = get_nhood_3D(shell_radius, [], thickness);   
    
     buf2 = zeros(movie_size); 
    for kkk = 1 : length(r1_shift)
        buf2(max(min(center_A_in(1)+ r1_shift(kkk), movie_size(1)), 1), ... % slightly larger
             max(min(center_A_in(2)+ r2_shift(kkk), movie_size(2)), 1), ...
             max(min(center_A_in(3)+ r3_shift(kkk), movie_size(3)), 1)) = 1;               
    end 
end