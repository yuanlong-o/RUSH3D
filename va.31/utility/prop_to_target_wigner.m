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