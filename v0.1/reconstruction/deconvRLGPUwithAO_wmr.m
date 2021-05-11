function Xguess = deconvRLGPUwithAO_wmr(maxIter, Xguess, blur_image, psf, Nnum, f_shift, savepath, AOstar, Angle_range, Block_num, defocus, Nbx,Nby,num_block,threshhold,sidelobe)
%% Size
[img_r,img_c,~] = size(blur_image);
[psf_r,~,~,~,allz] = size(psf);
%% Axial parts
if img_r>=700 && num_block ~= 1
    index_c = round(linspace(1,allz+1,num_block));
else
    index_c = [1,allz+1];
end
%% Weight for Each Angle
weight = squeeze(sum(sum(sum(blur_image,1),2),5));
for i = 1:Nnum
    for j = 1:Nnum
        if  ((i-8)^2+(j-8)^2 > Angle_range)
            weight(i,j) = 0;
        end
    end
end
weight=squeeze(weight./sum(weight(:)));
iter_weight = 0.8*1/max(weight(:)); %%最大迭代权重
%%
load('./seq15.mat');
blur_image = single(blur_image);
weight1 = gpuArray(single(weight));

for i=1:maxIter
    if AOstar == 1
        
        map_wavshape=zeros(Nnum,Nnum,Nbx,Nby,2);
        map_wavshape_weight=ones(Nnum,Nnum,Nbx,Nby);
        
        N3=round( 0.85*size(blur_image,1)/(Nbx)/2 )*2+1;
        borderx=(size(blur_image,1)-N3*Nbx)/2;
        bordery=(size(blur_image,2)-N3*Nby)/2;
        [coordinate1,coordinate2]=meshgrid(1:size(blur_image,1),1:size(blur_image,2));
        x=[-7:7];
        [Sx,Sy]=meshgrid(x,x);
        for u=1:Nnum
            for v=1:Nnum
                if ((u-8)^2+(v-8)^2)>Angle_range
                    Sx(u,v)=0;
                    Sy(u,v)=0;
                end
            end
        end
        mask_wjm=zeros(Nnum,Nnum);
        [xx,yy]=meshgrid(x,x);
        mask_wjm(xx.^2+yy.^2<=Angle_range)=1;
        
        if i>1
            for u_2=1:Nnum
                for v_2=1:Nnum
                    [v,u] = find((u_2-1)*Nnum+v_2 == seq);
                    if weight1(u,v)==0 || mask_wjm(u,v)==0
                        continue;
                    else
                        %%middle-Zslice FP
                        sumupXG1 = gpuArray.zeros(img_r,img_c,'single');
                        for block = max(round(0.2*length(index_c)),1):round(0.8*length(index_c))-1
                            cstart = index_c(block);
                            cend = index_c(block+1)-1;
                            tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                            tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                            tmp2 = fft3gyd(tmp2);
                            sumupXG1 = sumupXG1+sum(fft3gyd(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                        end
                        clear tmp2;
                        sumupXG1=abs(ifftshift(ifftn(sumupXG1)))./allz;
                        %%side-Zslice FP
                        sumupXG2 = gpuArray.zeros(img_r,img_c,'single');
                        for block = [1:max(round(0.2*length(index_c)),1)-1,round(0.8*length(index_c)):length(index_c)-1]
                            cstart = index_c(block);
                            cend = index_c(block+1)-1;
                            tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                            tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                            tmp2 = fft3gyd(tmp2);
                            sumupXG2 = sumupXG2+sum(fft3gyd(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                        end
                        clear tmp2;
                        sumupXG2=abs(ifftshift(ifftn(sumupXG2)))./allz;
                        for uu=1:Nbx
                            for vv=1:Nby
                                
                                sub_HXguess1=sumupXG1(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                                sub_HXguess2=sumupXG2(borderx+(uu-1)*N3+1:borderx+uu*N3,bordery+(vv-1)*N3+1:bordery+vv*N3);
                                sub_blur_image=gpuArray(squeeze(blur_image(borderx+(uu-1)*N3+1-sidelobe:borderx+uu*N3+sidelobe,bordery+(vv-1)*N3+1-sidelobe:bordery+vv*N3+sidelobe,u,v)));
                                corr_map=gather(normxcorr2(sub_HXguess1+sub_HXguess2,sub_blur_image));
                                [testa,testb]=find(corr_map==max(corr_map(:)));
                                %                                 f1 = fit((1:2*N3+2*sidelobe-1)',corr_map(testa(1),:)','gauss1','Lower',[0,testb(1)-1,-inf],'Upper',[1,testb(1)+1,inf]);
                                %                                 f2 = fit((1:2*N3+2*sidelobe-1)',corr_map(:,testb(1)),'gauss1','Lower',[0,testa(1)-1,-inf],'Upper',[1,testa(1)+1,inf]);
                                %                                 testa=f2.b1;
                                %                                 testb=f1.b1;
                                map_wavshape(u,v,uu,vv,1)=testa-size(sub_blur_image,1)+sidelobe;
                                map_wavshape(u,v,uu,vv,2)=testb-size(sub_blur_image,2)+sidelobe; %%似乎有点问题
                                disp(['uu = ',num2str(uu),'/',num2str(Nbx),',vv = ',num2str(vv),'/',num2str(Nby),' has been AOed!']);
                            end
                        end
                       clear sumupXG1 sumupXG2 sub_blur_image sub_HXguess1 sub_HXguess2
                    end
                end
                
            end
            %%去散焦和倾斜
            for uu=1:Nbx
                for vv=1:Nby
                    cx=map_wavshape(8,8,uu,vv,1);cy=map_wavshape(8,8,uu,vv,2);
                    map_wavshape(:,:,uu,vv,1)=(squeeze(map_wavshape(:,:,uu,vv,1))-cx).*mask_wjm;
                    map_wavshape(:,:,uu,vv,2)=(squeeze(map_wavshape(:,:,uu,vv,2))-cy).*mask_wjm;
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
                        k1 = Sy.*squeeze(map_wavshape(:,:,uu,vv,1)).*mask_wjm+Sx.*squeeze(map_wavshape(:,:,uu,vv,2)).*mask_wjm;
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
            save(strcat(savepath,'map_wavshape_test_scan_iter',num2str(i),'_block_',num2str(Block_num),'.mat'),'map_wavshape');            
        end    
    end
    
    for u_2=1:Nnum
        for v_2=1:Nnum
            [v,u] = find((u_2-1)*Nnum+v_2 == seq); % Iterate from the outside to the inside
            
            if weight1(u,v)==0
                continue;
            else
                if AOstar>0
                  %% AO分块连续
                    map_wavshape_x=squeeze(map_wavshape(u,v,:,:,1));
                    map_wavshape_y=squeeze(map_wavshape(u,v,:,:,2));
                    map_wavshape_xx=imresize(map_wavshape_x,[size(blur_image,1),size(blur_image,2)],'bilinear');
                    map_wavshape_yy=imresize(map_wavshape_y,[size(blur_image,1),size(blur_image,2)],'bilinear');
                    blur_image_uv=gpuArray(interp2(coordinate1,coordinate2,blur_image(:,:,u,v),coordinate1+map_wavshape_yy,coordinate2+map_wavshape_xx,'cubic',0));
                else
                    blur_image_uv=gpuArray(blur_image(:,:,u,v));
                end
                sumupXG = gpuArray.zeros(img_r,img_c,'single');
                for block = 1:length(index_c)-1
                    cstart = index_c(block);
                    cend = index_c(block+1)-1;
                    tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                    tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = flip(gpuArray(squeeze(psf(:,:,u,v,cstart:cend))),3);
                    tmp2 = fft3gyd(tmp2);
                    sumupXG = sumupXG+sum(fft3gyd(gpuArray(single(Xguess(:,:,cstart:cend)))).*tmp2,3);
                end
                clear tmp2;
                sumupXG=sumupXG./allz;
                for block = 1:length(index_c)-1
                    cstart = index_c(block);
                    cend = index_c(block+1)-1;
                    tmp2 = gpuArray.zeros(img_r,img_c,cend-cstart+1,'single');
                    tmp2((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,:)  = gpuArray(squeeze(psf(:,:,u,v,cstart:cend)));
                    tmp2 = fft3gyd(rot90(tmp2,2));
                    HXguessBack = (ifft3gyd(repmat(sumupXG,[1,1,cend-cstart+1]).*tmp2));
                    errorBack = ifft3gyd(fftshift(fftn(ifftshift(blur_image_uv))).*tmp2);
                    clear tmp2;
                    errorBack = real(errorBack./HXguessBack);
                    Xguess(:,:,cstart:cend) = gather((gpuArray(single(Xguess(:,:,cstart:cend))).*errorBack*weight1(u,v)*iter_weight+...
                        (1-weight1(u,v)*iter_weight).*gpuArray(single(Xguess(:,:,cstart:cend)))));
                end
                clear blur_image_uv errorBack HXguessBack sumupXG; 
                disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
            end
        end
    end
    A = Xguess(sidelobe+1:end-sidelobe,sidelobe+1:end-sidelobe,:);
    if size(A,1) > 2000 && size(A,3) >100
        imwriteTFSK(single(A(:,:,1:round(size(A,3)/2))), ...
            [savepath,'AO_',num2str(AOstar),'_iter_',num2str(i),'.0.tiff']);
        imwriteTFSK(single(A(:,:,round(size(A,3)/2)+1:end)), ...
            [savepath,'AO_',num2str(AOstar),'_iter_',num2str(i),'.1.tiff']);
    else
        imwriteTFSK(single(A), [savepath,'AO_',num2str(AOstar),'_iter_',num2str(i),'.tiff']);
    end
end
if(f_shift)
    Xguess = fftshift(Xguess,1);
    Xguess = fftshift(Xguess,2);
    Xguess=gather(Xguess);
else
    Xguess=gather(Xguess);
end
end

function output = fft3gyd(input)
output = fftshift(fftn(ifftshift(input)));
end

function output = ifft3gyd(input)
output = fftshift(ifftn(ifftshift(input)));
end
