function [Xguess] = recon_module(recon_param, Xguess, psf, WDF, ...
                            recon_name_perfix ,view_array, frame)
%% Reconstruction Module reconstuction with DAO
% This program is used to reconstruct volume from wdf with digital
% AOstar = 2
% recon_angle_range
% abberation optics correction
% Last update: 05/16/2021. MW
% Last update: 03/18/2022. MW

% Input:
% psf_param          Nnum
% recon_param       
% Xguess             Initialized X volume
% psf                N x N x V x Z
% WDF                M x M x V  
% recon_savepath
% recon_name_perfix
% frame              use to name the reconstruction file

% Output:
% Xguess
% output file        map_waveshape and Xguess


Nnum = recon_param.Nnum;
writemode = recon_param.writemode;
dispmode = recon_param.dispmode;
view_range = recon_param.view_range; % About 25 Angle views within the iteration
AOstar = recon_param.AOstar; % 1 for DAO; 0 for no DAO 
maxIter = recon_param.maxIter; % Max iteration times: 2 or 3 is enough
defocus = recon_param.defocus; % 1 for defocus, 0 for no defocus
Nbx = recon_param.Nbx; 
Nby = recon_param.Nby; % Block apart 5 x 5 pieces when DAO 
threshhold = recon_param.threshhold; % Shift should not be allowed to exceed [-25,25]
margin = recon_param.margin; % margin overlap
% new parameter
recon_angle_range = recon_param.angle_range;
if AOstar == 2
    restart = recon_param.restart;
else
    restart = 20;
end


[img_r,img_c,~] = size(WDF);
[psf_r,~,~,psf_z] = size(psf);
view_num = length(view_array);
weight = squeeze(sum(sum(sum(WDF,1),2),4));
weight=squeeze(weight./sum(weight(:)));
iter_weight = 0.4*1/max(weight(:)); % iteration weight
cen_num = ceil(Nnum/2); 

[~,~,seq_ind] = Spiral_circle(Nnum,view_range);
WDF = single(WDF);
weight1 = gpuArray(single(weight));
Xguess=gpuArray(single(Xguess));
if AOstar == 2
    Xguess_tmp = Xguess;
end
%% start iteration
for i=1:maxIter
    
    if AOstar > 0
        map_wavshape=zeros(Nnum,Nnum,Nbx,Nby,2);
        map_wavshape_weight=ones(Nnum,Nnum,Nbx,Nby);
        
        N3=round( 0.88*size(WDF,1)/(Nbx)/2 )*2+1;
        borderx=round((size(WDF,1)-N3*Nbx)/2);
        bordery=round((size(WDF,2)-N3*Nby)/2);
        [coordinate1,coordinate2]=meshgrid(1:size(WDF,2),1:size(WDF,1));
        x = -floor(Nnum/2):floor(Nnum/2);
        [Sx,Sy]=meshgrid(x,x);
        mask = Sx.^2 + Sy.^2 <= view_range^2;
        Sx = Sx.*mask;
        Sy = Sy.*mask;
        
        if i>1
            for v_ind=1:view_num
                view = seq_ind(v_ind);
                u = view_array{view}(1);
                v = view_array{view}(2);
                tmp2 = gpuArray.zeros(img_r,img_c,psf_z,'single');
                tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,view ,:))),3);
                tmp2 = fft3_new(tmp2);
                sumupXG1 = single(sum(fft3_new(Xguess).*tmp2,3));
                sumupXG1=abs(ifftshift(ifftn(sumupXG1)))./psf_z;
                
                clear tmp2;
                for uu=1:Nbx
                    for vv=1:Nby
                        
                        sub_HXguess=sumupXG1(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                        sub_blur_image=gpuArray(squeeze(WDF(borderx+(uu-1)*N3+1-margin:borderx+uu*N3+margin,bordery+(vv-1)*N3+1-margin:bordery+vv*N3+margin,view)));
                        corr_map=gather(normxcorr2(sub_HXguess,sub_blur_image));
                        [shift_a,shift_b]=find(corr_map==max(corr_map(:)));
                        
                        map_wavshape(u,v,uu,vv,1)=shift_a(1)-size(sub_blur_image,1)+margin;
                        map_wavshape(u,v,uu,vv,2)=shift_b(1)-size(sub_blur_image,2)+margin;
                        if dispmode == 1
                            disp(['uu = ',num2str(uu),'/',num2str(Nbx),',vv = ',num2str(vv),'/',num2str(Nby),' has been AOed!']);
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
                            map_wavshape(u,v,uu,vv,1)=min(max(map_wavshape(u,v,uu,vv,1).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);
                            map_wavshape(u,v,uu,vv,2)=min(max(map_wavshape(u,v,uu,vv,2).*map_wavshape_weight(u,v,uu,vv),-threshhold),threshhold);
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
    
    if i == restart && AOstar == 2
        Xguess = Xguess_tmp;
    end
    
    for v_ind = 1: view_num
        tic;
        view = seq_ind(v_ind);
        u = view_array{view}(1);
        v = view_array{view}(2);
        if (u-cen_num)^2 + (v-cen_num)^2 >= recon_angle_range
            continue
        else
            if AOstar>0
                % AO correction
                map_wavshape_x=squeeze(map_wavshape(u,v,:,:,1));
                map_wavshape_y=squeeze(map_wavshape(u,v,:,:,2));
                map_wavshape_xx=imresize(map_wavshape_x,[size(WDF,1),size(WDF,2)],'bilinear');
                map_wavshape_yy=imresize(map_wavshape_y,[size(WDF,1),size(WDF,2)],'bilinear');
                blur_image_uv=gpuArray(interp2(coordinate1,coordinate2,WDF(:,:,view),coordinate1+map_wavshape_yy,coordinate2+map_wavshape_xx,'cubic',0));
            else
                blur_image_uv=gpuArray(WDF(:,:,view));
            end
            
            tmp2 = gpuArray.zeros(img_r,img_c,psf_z,'single');
            tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:)  = gpuArray(squeeze(psf(:,:,view,:)));
            tmp3 = flip(tmp2,3);
            tmp3 = fft3_new(tmp3);
            sumupXG = single(sum(fft3_new(Xguess).*tmp3,3));
            sumupXG = sumupXG./psf_z;
            
            tmp2 = fft3_new(rot90(tmp2,2));
            HXguessBack = ifft3_new(repmat(sumupXG,[1,1,psf_z]).*tmp2);
            errorBack = ifft3_new(fftshift(fftn(ifftshift(blur_image_uv))).*tmp2);
            errorBack = real(errorBack./HXguessBack);
            clear HXguessBack tmp2 tmp3;
            Xguess = Xguess.*errorBack*weight1(view)*iter_weight+(1-weight1(view)*iter_weight).*Xguess;
            clear errorBack sumupXG;
            tt = toc;
            if dispmode == 1
                disp(['iter ' num2str(i) ' | ' num2str(maxIter) ', (v=' num2str(view)  '),  Energy=' num2str(sum(Xguess(:))) ', takes' num2str(tt) 's']);
            end
        end
    end
    
    if writemode == 1 || (writemode == 2 && i == maxIter)
        A = abs(double(gather(Xguess)));
        imwriteTFSK(single(A) , [recon_name_perfix,'_vid', num2str(frame),'_iter_',num2str(i),'.tiff']);
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