function [Xguess] = reconstruction_module(psf_param, recon_param, Xguess, psf, WDF, ...
                            recon_savepath, recon_name_perfix, frame)
%% Reconstruction Module reconstuction with DAO
% This program is used to reconstruct volume from wdf with digital
% abberation optics correction
% Last update: 05/16/2021. MW

% Input:
% psf_param          Nnum
% recon_param       
% Xguess             Initialized X volume
% psf                N x N x Nnum x Nnum x Z
% WDF                M x M x Nnum x Nnum  
% recon_savepath
% recon_name_perfix
% frame              use to name the reconstruction file

% Output:
% Xguess
% output file        map_waveshape and Xguess


Nnum = psf_param.Nnum;
angle_range = recon_param.angle_range; % About 25 Angle views within the iteration
AOstar = recon_param.AOstar; % 1 for DAO; 0 for no DAO 
maxIter = recon_param.maxIter; % Max iteration times: 2 or 3 is enough
defocus = recon_param.defocus; % 1 for defocus, 0 for no defocus
Nbx = recon_param.Nbx; 
Nby = recon_param.Nby; % Block apart 5 x 5 pieces when DAO 
num_block = recon_param.num_block; % Axial block for 10 when forword propagate 
threshhold = recon_param.threshhold; % Shift should not be allowed to exceed [-25,25]
margin = recon_param.margin; % margin overlap

[img_r,img_c,~] = size(WDF);
[psf_r,~,~,~,allz] = size(psf);

if img_r>=700 && num_block ~= 1
    index_c = round(linspace(1,allz+1,num_block));
else
    index_c = [1,allz+1];
end

weight = squeeze(sum(sum(sum(WDF,1),2),5));
for i = 1:Nnum
    for j = 1:Nnum
        if  ((i-8)^2+(j-8)^2 > angle_range)
            weight(i,j) = 0;
        end
    end
end
weight=squeeze(weight./sum(weight(:)));
iter_weight = 0.8*1/max(weight(:)); % iteration weight

load('./reconstruction_module/seq15.mat');
WDF = single(WDF);
weight1 = gpuArray(single(weight));

%% start iteration
for i=1:maxIter
    if AOstar == 1
        map_wavshape=zeros(Nnum,Nnum,Nbx,Nby,2);
        map_wavshape_weight=ones(Nnum,Nnum,Nbx,Nby);
        
        N3=round( 0.85*size(WDF,1)/(Nbx)/2 )*2+1;
        borderx=(size(WDF,1)-N3*Nbx)/2;
        bordery=(size(WDF,2)-N3*Nby)/2;
        [coordinate1,coordinate2]=meshgrid(1:size(WDF,1),1:size(WDF,2));
        x = -7:7;
        [Sx,Sy]=meshgrid(x,x);
        for u=1:Nnum
            for v=1:Nnum
                if ((u-8)^2+(v-8)^2)>angle_range
                    Sx(u,v)=0;
                    Sy(u,v)=0;
                end
            end
        end
        mask=zeros(Nnum,Nnum);
        [xx,yy]=meshgrid(x,x);
        mask(xx.^2+yy.^2<=angle_range)=1;
        
        
        if i>1
            for u_2=1:Nnum
                for v_2=1:Nnum
                    [v,u] = find((u_2-1)*Nnum+v_2 == seq);
                    if weight1(u,v)==0 || mask(u,v)==0
                        continue;
                    
                    % DAO calculate correlation
                    else
                        sumupXG1 = gpuArray.zeros(img_r,img_c,'single');
                        for block = max(round(0.2*length(index_c)),1):round(0.8*length(index_c))-1
                            cstart = index_c(block);
                            cend = index_c(block+1)-1;
                            tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                            tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                            tmp2 = fft3_new(tmp2);
                            sumupXG1 = sumupXG1+sum(fft3_new(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                        end
                        sumupXG1=abs(ifftshift(ifftn(sumupXG1)))./allz;
                        sumupXG2 = gpuArray.zeros(img_r,img_c,'single');
                        clear tmp2;
                        for block = [1:max(round(0.2*length(index_c)),1)-1,round(0.8*length(index_c)):length(index_c)-1]
                            cstart = index_c(block);
                            cend = index_c(block+1)-1;
                            tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                            tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                            tmp2 = fft3_new(tmp2);
                            sumupXG2 = sumupXG2+sum(fft3_new(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                        end
                        sumupXG2=abs(ifftshift(ifftn(sumupXG2)))./allz;
                        clear tmp2;
                        for uu=1:Nbx
                            for vv=1:Nby
                                
                                sub_HXguess1=sumupXG1(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                                sub_HXguess2=sumupXG2(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                                sub_blur_image=gpuArray(squeeze(WDF(borderx+(uu-1)*N3+1-margin:borderx+uu*N3+margin,bordery+(vv-1)*N3+1-margin:bordery+vv*N3+margin,u,v)));
                                corr_map=gather(normxcorr2(sub_HXguess1+sub_HXguess2,sub_blur_image));
                                [shift_a,shift_b]=find(corr_map==max(corr_map(:)));

                                map_wavshape(u,v,uu,vv,1)=shift_a-size(sub_blur_image,1)+margin;
                                map_wavshape(u,v,uu,vv,2)=shift_b-size(sub_blur_image,2)+margin;
                                disp(['uu = ',num2str(uu),'/',num2str(Nbx),',vv = ',num2str(vv),'/',num2str(Nby),' has been AOed!']);
                            end
                        end
                    end
                end
            end
            clear sub_HXguess1 sub_HXguess2 sub_blur_image sumupXG1 sumupXG2
            
            for uu=1:Nbx
                for vv=1:Nby
                    cx=map_wavshape(8,8,uu,vv,1);cy=map_wavshape(8,8,uu,vv,2);
                    map_wavshape(:,:,uu,vv,1)=(squeeze(map_wavshape(:,:,uu,vv,1))-cx).*mask;
                    map_wavshape(:,:,uu,vv,2)=(squeeze(map_wavshape(:,:,uu,vv,2))-cy).*mask;
                end
            end
            for uu=1:Nbx
                for vv=1:Nby
                    for u=1:Nnum
                        for v=1:Nnum
                            map_wavshape(u,v,uu,vv,1)=min(max(map_wavshape(u,v,uu,vv,1).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);%+min(max(cx,-10),10);
                            map_wavshape(u,v,uu,vv,2)=min(max(map_wavshape(u,v,uu,vv,2).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);%+min(max(cx,-10),10);
                        end
                    end
                end
            end
            if defocus == 1
                for uu=1:Nbx
                    for vv=1:Nby
                        k1 = Sy.*squeeze(map_wavshape(:,:,uu,vv,1)).*mask+Sx.*squeeze(map_wavshape(:,:,uu,vv,2)).*mask;
                        k2 = Sx.*Sx+Sy.*Sy;
                        k=sum(k1(:))/sum(k2(:));
                        map_wavshape(:,:,uu,vv,1)=squeeze(map_wavshape(:,:,uu,vv,1))-k*Sy;
                        map_wavshape(:,:,uu,vv,2)=squeeze(map_wavshape(:,:,uu,vv,2))-k*Sx;
                        
                        for u=1:Nnum
                            for v=1:Nnum
                                map_wavshape(u,v,uu,vv,1)=min(max(map_wavshape(u,v,uu,vv,1).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);%+min(max(cx,-10),10);
                                map_wavshape(u,v,uu,vv,2)=min(max(map_wavshape(u,v,uu,vv,2).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);%+min(max(cx,-10),10);
                            end
                        end
                        
                    end
                end
            end
            save(strcat(recon_name_perfix,'map_wavshape_iter',num2str(i),'.mat'),'map_wavshape');            
        end    
    end
    
    for u_2=1:Nnum
        for v_2=1:Nnum
            [v,u] = find((u_2-1)*Nnum+v_2 == seq);
            
            if weight1(u,v)==0
                continue;
            else
                if AOstar>0
                 % AO correction
                    map_wavshape_x=squeeze(map_wavshape(u,v,:,:,1));
                    map_wavshape_y=squeeze(map_wavshape(u,v,:,:,2));
                    map_wavshape_xx=imresize(map_wavshape_x,[size(WDF,1),size(WDF,2)],'bilinear');
                    map_wavshape_yy=imresize(map_wavshape_y,[size(WDF,1),size(WDF,2)],'bilinear');
                    blur_image_uv=gpuArray(interp2(coordinate1,coordinate2,WDF(:,:,u,v),coordinate1+map_wavshape_yy,coordinate2+map_wavshape_xx,'cubic',0));
                else
                    blur_image_uv=gpuArray(WDF(:,:,u,v));
                end
                sumupXG = gpuArray.zeros(img_r,img_c,'single');
                for block = 1:length(index_c)-1
                    cstart = index_c(block);
                    cend = index_c(block+1)-1;
                    tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                    tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                    tmp2 = fft3_new(tmp2);
                    sumupXG = sumupXG+sum(fft3_new(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                end
                sumupXG=sumupXG./allz;
                for block = 1:length(index_c)-1
                    cstart = index_c(block);
                    cend = index_c(block+1)-1;
                    tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                    tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = gpuArray(squeeze(psf(:,:,u,v,cstart:cend)));
                    tmp2 = fft3_new(rot90(tmp2,2));
                    HXguessBack = (ifft3_new(repmat(sumupXG,[1,1,cend-cstart+1]).*tmp2));
                    errorBack = ifft3_new(fftshift(fftn(ifftshift(blur_image_uv))).*tmp2);
                    errorBack = real(errorBack./HXguessBack);
                    clear HXguessBack tmp2;
                    Xguess(:,:,cstart:cend) = gather((gpuArray(single(Xguess(:,:,cstart:cend))).*errorBack*weight1(u,v)*iter_weight+...
                        (1-weight1(u,v)*iter_weight).*gpuArray(single(Xguess(:,:,cstart:cend)))));
                    clear errorBack;
                end
                clear sumupXG;
                disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
            end
        end
    end
    A = Xguess(margin+1:end-margin,margin+1:end-margin,:);
    
    % save (if it is too large, save as several stacks)
    if size(A,1) > 2000 && size(A,3) > 100
        saveastiff(im2uint16(single(A(:,:,1:round(size(A,3)/2))) / max(A(:))), ...
            [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.0.tiff']);
        saveastiff(im2uint16(single(A(:,:,round(size(A,3)/2)+1:end)) / max(A(:))), ...
            [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.1.tiff']);
    else
        saveastiff(im2uint16(single(A) / max(A(:))), [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.tiff']);
    end
end

Xguess=gather(Xguess);

end

function output = fft3_new(input)
output = fftshift(fftn(ifftshift(input)));
end

function output = ifft3_new(input)
output = fftshift(ifftn(ifftshift(input)));
end
