function projection_wdf = forwardProjectAllViewConv(psf_param,psf,Xguess)


Nnum = psf_param.Nnum;

% Xguess=gpuArray(single(Xguess));
% psf = gpuArray(single(psf));

[~,~,~,~,psf_z] = size(psf);
[img_r,img_c,~] = size(Xguess);



projection_wdf = gpuArray.zeros(img_r,img_c,Nnum,Nnum,'single');

for u = 1:Nnum
    parfor v = 1:Nnum
        projection = single(zeros(img_r,img_c));
        for z = 1: psf_z
            projection = projection + conv2(Xguess(:,:,z),psf(:,:,u,v,z),'same');      
        end
        projection_wdf(:,:,u,v) = projection;
        disp(['  calc WDF (u=',num2str(u), ',  v=',num2str(v)  ]);
    end
end
projection_wdf = gather(projection_wdf);
end




