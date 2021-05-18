%%本文件已废弃 放弃更新 请勿使用



tic
startFrame = 18; %起始帧号(0~N-1)
frameJump = 1; %跳帧：针对多色和希望跳帧的情况（1~N)
name = 'mouse_liver_8_9hz_3x3_100.0ms_Full_Hardware_LaserCount1_190903213802.0.tiff';%输入的第一组tiff的名字（会自动往下找，请勿重命名原文件名）
realignName = 'realign/test';%输出图像前缀（请注意自己新建文件夹，否则没有输出，不会报错）
mkdir('realign');%新建输出的文件夹
Nshift = 3;%拍摄的扫描
Nx = 75;%取的微透镜阵列数 Nx = (pick - 1) /2
Ny = 75;%同上，Y方向，代码有问题，建议Nx==Ny
Nnum = 13;%微透镜阵列后的pixel数，一般不需要动
resize = 0;%是否resize到13x13输出
confName = './3x3.conf.sk.png';%扫描配置文件，需要根据Nshift做修改
groupCount = 40;%计算多少组
groupMode = 0;%0=间隔模式(1-9,10-18,...);1=连续模式(1-9,2-10,...)
autoCenterMode = 1;%自动找中心点
autoCenterFrame = startFrame;%自动找中心点用的那张图，默认取起始图片(zgx mode必须取扫描第一张），从0到N-1
centerX=1131;%中心点X，自动模式下失效，从0到N-1（在c++坐标系，从0开始！）
centerY=1087;%中心点Y，自动模式下失效，从0到N-1（在c++坐标系，从0开始！）
realignMode = 'LZ';%realignMode = 'ZGX'; %LZ mode 或者 ZGX mode
centerView = 'realign/centerView.only.tiff';%输出中心视角模式，若不需要，此处填一个'None'，需要则填tiff的名字， 若.tiff前有个.only，则只会输出中心视角（例如：outdir/a.only.tiff）
rotation =  0;%旋转光场图，顺时针，可填：0,90,180,270


if(autoCenterMode == 1)
    img = double(imread(name, autoCenterFrame + 1));
    kernal = fspecial('gaussian',[Nnum,Nnum],3);
    fed = imfilter(img, kernal);
    locMatrix = zeros(Nnum, Nnum);
    for i = 1:Nnum
        for j = 1:Nnum
            picMat = fed(i:Nnum:end, j:Nnum:end);
            locMatrix(i,j) = mean(mean(picMat));
        end
    end
    [yc,xc] = find(locMatrix == max(max(locMatrix)));
    centerXX = xc + floor(size(img,2) / Nnum / 2) * Nnum - 1;
    centerYY = yc + floor(size(img,1) / Nnum / 2) * Nnum - 1;
    fprintf('AutoCenter found center point x = %d, y = %d (range from 0~[size-1]), check at ImageJ\n', centerXX, centerYY);
    if(mod(centerXX, Nnum) ~= mod(centerX, Nnum) || mod(centerYY, Nnum) ~= mod(centerY, Nnum))
        warning('AutoCenter Point not match! Found(modded) = (%d, %d), input = (%d, %d), both at c++ range',...
            mod(centerXX, Nnum), mod(centerYY, Nnum), mod(centerX, Nnum), mod(centerY, Nnum));
    end
    locMatrixSort = sort(reshape(locMatrix,[1,Nnum * Nnum]));
    yyt = 0; xxt = 0;
    for i = 1:4
        [yyyt, xxxt] = find(locMatrix == locMatrixSort(i));
        yyt = yyyt / 4 + yyt;
        xxt = xxxt / 4 + xxt;
    end
    %TODO:Here we are facing a problem, when 2 least pixel reach the edge!
    if(abs(abs(yyt - yc) - Nnum / 2) > 1e-3 || abs(abs(xxt - xc) - Nnum / 2) > 1e-3 || min(min(locMatrix)) / max(max(locMatrix)) > 0.9)
        warning('呜呜呜，我找的中心点可能有问题……你能不能再检查一下用于找中心点的图呢QAQ');
    else
        ctemp = (0.9 - min(min(locMatrix)) / max(max(locMatrix))) / 0.2;
        if(ctemp > 1) ctemp = 1; end
        fprintf('AutoCenter calculation credibility = %.2f%%\n', ctemp * 100); %这个概率我瞎写的，大概就是最小值/最大值=0.9->0%（信噪比过低）;=0.7->100%（信噪比足够）
    end
    centerX = centerXX;
    centerY = centerYY;
    %clear xc xxt xxxt yc yyt yyyt locMatrixSort img kernal fed ctemp
    %t = imresize(locMatrix, [1300,1300], 'nearest');imshow(t,[]);
end



%old version
% if(autoCenterMode == 1)
%     img = imread(name, startFrame + 1);
%     sumX = sum(img);
%     sumY = sum(img,2);
%     sumX = sumX + circshift(sumX, 1, 2);
%     sumY = sumY + circshift(sumY, 1, 1);
%     centerX = round(size(img, 2) / 2 - Nnum / 2);
%     for i = round(size(img, 2) / 2 - Nnum / 2) : round(size(img, 2) / 2 - Nnum / 2) + Nnum
%         if(sumX(i) < sumX(centerX))
%             centerX = i;
%         end
%     end
%     centerY = round(size(img, 1) / 2 - Nnum / 2);
%     for i = round(size(img, 1) / 2 - Nnum / 2) : round(size(img, 1) / 2 - Nnum / 2) + Nnum
%         if(sumY(i) < sumY(centerY))
%             centerY = i;
%         end
%     end
%     centerX = centerX + floor(Nnum / 2) - 1;
%     centerY = centerY + floor(Nnum / 2) - 1;
%     fprintf('found center point x = %d, y = %d, check at ImageJ', centerX, centerY);
% end



%%
% ./ReAlign [CaptureNshift] [OutputFile] [RawData]
% [RawDataIdx] [CenterPixelX] [CenterPixelY]
% [PickLensSizeX] (DataGroup = 1) (GroupMode = 0){0=individual, 1=continues} (ShiftConfigFile = Auto)
% (FrameJump = 1) (Resize = 1) (ResizeShift = 13) (PickLensSizeY = X)
command = sprintf('./ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d %s %s %d',...
Nshift, realignName, name, startFrame, centerX, centerY, Nx,...
groupCount, groupMode, confName, frameJump, resize, Nnum, Ny, realignMode, centerView, rotation);
system(command);
toc