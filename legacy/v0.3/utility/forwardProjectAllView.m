function projection = forwardProjectAllView(psf_param,psf,Xguess)
%% This program is used for calculate forward projection
% input
% psf             M x M x z
% Xguess          N x N x z
% output
% projection      M x M

Nnum = psf_param.Nnum;

Xguess=gpuArray(single(Xguess));
psf = gpuArray(single(psf));

[psf_r,~,~,~,psf_z] = size(psf);
[img_r,img_c,~] = size(Xguess);



projection = gpuArray.zeros(img_r,img_c,Nnum,Nnum,'single');
for u = 1:Nnum
    for v = 1:Nnum
        tmp2 = gpuArray.zeros(img_r,img_c,psf_z,'single');
        tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:)  = flip(squeeze(psf(:,:,u,v,:)),3);
        tmp2 = fft3_new(tmp2);
        sumupXG1 = single(sum(fft3_new(Xguess).*tmp2,3));
        clear tmp2;
        projection(:,:,u,v)=abs(ifftshift(ifftn(sumupXG1)))./psf_z; 
        disp(['  calc WDF (u=',num2str(u), ',  v=',num2str(v)  ]);
    end
end
projection = gather(projection);
end


function output = fft3_new(input)
output = fftshift(fftn(ifftshift(input)));
end
