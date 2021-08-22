function [Xguess] = reconstruction_module(recon_param, Xguess, psf, WDF, ...
                            recon_name_perfix, frame)
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


Nnum = recon_param.Nnum;
writemode = recon_param.writemode;
dispmode = recon_param.dispmode;
angle_range = recon_param.angle_range; % About 25 Angle views within the iteration
AOstar = recon_param.AOstar; % 1 for DAO; 0 for no DAO 
maxIter = recon_param.maxIter; % Max iteration times: 2 or 3 is enough
defocus = recon_param.defocus; % 1 for defocus, 0 for no defocus
Nbx = recon_param.Nbx; 
Nby = recon_param.Nby; % Block apart 5 x 5 pieces when DAO 
threshhold = recon_param.threshhold; % Shift should not be allowed to exceed [-25,25]
margin = recon_param.margin; % margin overlap

[img_r,img_c,~] = size(WDF);
[psf_r,~,~,~,psf_z] = size(psf);

weight = squeeze(sum(sum(sum(WDF,1),2),5));
for i = 1:Nnum
    for j = 1:Nnum
        if  ((i-ceil(Nnum/2))^2+(j-ceil(Nnum/2))^2)>angle_range
            weight(i,j) = 0;
        end
    end
end
weight=squeeze(weight./sum(weight(:)));
iter_weight = 0.8*1/max(weight(:)); % iteration weight


seq = double(Spiral(Nnum));
WDF = single(WDF);
weight1 = gpuArray(single(weight));
Xguess=gpuArray(single(Xguess));
%% start iteration
for i=1:maxIter
    
    if AOstar == 1
        map_wavshape=zeros(Nnum,Nnum,Nbx,Nby,2);
        map_wavshape_weight=ones(Nnum,Nnum,Nbx,Nby);
        
        N3=round( 0.85*size(WDF,1)/(Nbx)/2 )*2+1;
        borderx=round((size(WDF,1)-N3*Nbx)/2);
        bordery=round((size(WDF,2)-N3*Nby)/2);
        [coordinate1,coordinate2]=meshgrid(1:size(WDF,2),1:size(WDF,1));
        x = -floor(Nnum/2):floor(Nnum/2);
        [Sx,Sy]=meshgrid(x,x);
        for u=1:Nnum
            for v=1:Nnum
                if ((u-ceil(Nnum/2))^2+(v-ceil(Nnum/2))^2)>angle_range
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

                        tmp2 = gpuArray.zeros(img_r,img_c,psf_z,'single');
                        tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,:))),3);
                        tmp2 = fft3_new(tmp2);
                        sumupXG1 = single(sum(fft3_new(Xguess).*tmp2,3));
                        sumupXG1=abs(ifftshift(ifftn(sumupXG1)))./psf_z;
                        
                        clear tmp2;
                        for uu=1:Nbx
                            for vv=1:Nby
                                
                                sub_HXguess=sumupXG1(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                                sub_blur_image=gpuArray(squeeze(WDF(borderx+(uu-1)*N3+1-margin:borderx+uu*N3+margin,bordery+(vv-1)*N3+1-margin:bordery+vv*N3+margin,u,v)));
                                corr_map=gather(normxcorr2(sub_HXguess,sub_blur_image));
                                [shift_a,shift_b]=find(corr_map==max(corr_map(:)));

                                map_wavshape(u,v,uu,vv,1)=shift_a-size(sub_blur_image,1)+margin;
                                map_wavshape(u,v,uu,vv,2)=shift_b-size(sub_blur_image,2)+margin;
                                if dispmode == 1
                                    disp(['uu = ',num2str(uu),'/',num2str(Nbx),',vv = ',num2str(vv),'/',num2str(Nby),' has been AOed!']);
                                end
                            end
                        end
                    end
                end
            end
            clear sub_HXguess1 sub_HXguess2 sub_blur_image sumupXG1 sumupXG2
            
            for uu=1:Nbx
                for vv=1:Nby
                    cx=map_wavshape(ceil(Nnum/2),ceil(Nnum/2),uu,vv,1);cy=map_wavshape(ceil(Nnum/2),ceil(Nnum/2),uu,vv,2);
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
            if writemode == 1
                save(strcat(recon_name_perfix,'map_wavshape_iter',num2str(i),'.mat'),'map_wavshape');
            end
        end    
    end
    
    for u_2=1:Nnum
        for v_2=1:Nnum
            tic;
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


                tmp2 = gpuArray.zeros(img_r,img_c,psf_z,'single');
                tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:)  = gpuArray(squeeze(psf(:,:,u,v,:)));
                tmp3 = flip(tmp2,3);
                tmp3 = fft3_new(tmp3);
                
                sumupXG = single(sum(fft3_new(Xguess).*tmp3,3));
                
                sumupXG=sumupXG./psf_z;
                tmp2 = fft3_new(rot90(tmp2,2));
                HXguessBack = (ifft3_new(repmat(sumupXG,[1,1,psf_z]).*tmp2));
                errorBack = ifft3_new(fftshift(fftn(ifftshift(blur_image_uv))).*tmp2);
                errorBack = real(errorBack./HXguessBack);
                clear HXguessBack tmp2 tmp3;
                Xguess = Xguess.*errorBack*weight1(u,v)*iter_weight+(1-weight1(u,v)*iter_weight).*Xguess;
                clear errorBack;

                clear sumupXG;
                tt = toc;
                if dispmode == 1
                    disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
                    disp([num2str(tt),'s']);
                end
            end
        end
    end
    
    A = abs(double(gather(Xguess)));
    %A = A(margin+1:end-margin,margin+1:end-margin,:);
    
    if writemode == 1 || (writemode == 2 && i == maxIter)
        
        if size(A,1) > 2000 && size(A,3) > 100
            imwriteTFSK(single(A(:,:,1:round(size(A,3)/2))), ...
                [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.0.tiff']);
            imwriteTFSK(single(A(:,:,round(size(A,3)/2)+1:end)), ...
                [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.1.tiff']);
        else
            imwriteTFSK(single(A) , [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.tiff']);
        end
    end
end
%tt = toc;
Xguess=gather(Xguess);


end

function output = fft3_new(input)
output = fftshift(fftn(ifftshift(input)));
end

function output = ifft3_new(input)
output = fftshift(ifftn(ifftshift(input)));
end
