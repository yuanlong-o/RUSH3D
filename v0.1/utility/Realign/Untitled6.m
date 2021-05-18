clear;
fprintf('ver = 1.2.3\n');
% 更新(19-11-12)：
% 增加了LaserMergeOption
% 大致解决了图像数量统计不正确的问题(12001张图）
% 保证了当不输出中心视角时，自动Merge不会崩溃的问题
% 增加了try-catch语句以保证正确
% 修正了一个自动中心点处理过爆图像出错的问题
% 增加了中心点算法在边缘的可信度减少比例
% 更新（19-12-1）
% 修正了一个在指定中心点时会崩溃的bug
% 更新（20-1-2）
% 增加了preResize（实验阶段）
% 更新（20-3-28）
% 增加了手动FRAMEJUMP
% 修正了处理CV文件的bug
% 更新（20-5-22）
% 增加了SlightRotate


RawDir=fullfile('X:\zhu-ty\Realign');
Nnum = 13;%微透镜阵列后的pixel数，一般不需要动
resize = 0;%是否resize到13x13输出
groupMode = 1;%0=间隔模式(1-9,10-18,...);1=连续模式(1-9,2-10,...)
centerX=1015;%中心点X，自动模式下失效，从0到N-1（在c++坐标系，从0开始！） 请填preResize之前的中心点，程序中会自动乘上preResize
centerY=255;%中心点Y，自动模式下失效，从0到N-1（在c++坐标系，从0开始！）
realignMode = 'LZ';%realignMode = 'ZGX'; %LZ mode 或者 ZGX mode
rotation = 0;%旋转光场图，顺时针，可填：0,90,180,270
preResize = 2050/2048;%Resize光场图（预处理，默认值1.0）
slightRotation = 0;%轻微旋转光场图（插值） 这个轻微旋转在上面的旋转之前做，带插值
%顺序：先根据所填（找）的中心点做轻微旋转，再resize光场图，再resize中心点坐标，再裁剪图像，再进行大旋转，然后进入Realign抽取
%特别的，请注意，当rotation不等于0的时候，所搜寻的conf文件名应当符合下列格式
%NxN.conf.sk.rot[Angle].png  Sample:
%5x5.conf.sk.rot180.png
%13x13.conf.sk.rot180.png

FRAMEJUMP = 0;%=0 自动根据LaserCount决定分几个通道 =n 手动分通达

autoCenterMode = 1;%自动找中心点
centerviewOnly = 0;
LaserMergeOption = 1;%=1自动merge多通道，=0不merge
MergeContrastRatio = 2;




dirOutput=dir(fullfile(RawDir,'*.0.tiff'));
source={dirOutput.name}';

for sourceIdx = 1:size(source,1)
    if(source{sourceIdx}(end - 8) == 'C')
        continue;
    end
    try
        ticer = tic;
        front = source{sourceIdx}(1:end-6);
        dirOutput=dir(fullfile(RawDir,strcat(front, '*')));
        this_name = {dirOutput.name}';
        Itmp = imread(strcat(RawDir, '/',source{sourceIdx}), 1);
        [h,w] = size(Itmp);
        Nx = floor(w * preResize / (2*Nnum)) - 2;
        Ny = floor(h * preResize / (2*Nnum)) - 2;
        imgSize = h * w * 2;
        %此处改为分文件统计
        fileSizes = cell2mat({dirOutput.bytes});
        imgCount = 0;
        for fileIdx = 1:size(fileSizes, 2)
            imgCount = imgCount + floor(fileSizes(fileIdx) / imgSize);
        end
        %imgFullSize = sum(cell2mat({dirOutput.bytes}));
        %imgCount = floor(imgFullSize / imgSize);
        name = strcat(RawDir, '/',source{sourceIdx});

        %laserCountStrIdx = findstr(s,'wood')
        if(FRAMEJUMP == 0)
            [laserCountIdxS, laserCountIdxE] = regexp(front, '(?<=_LaserCount)[^/]+(?=_)');
            frameJump = str2num(front(laserCountIdxS:laserCountIdxE)); %跳帧：针对多色和希望跳帧的情况(1~N)
        else
            frameJump = FRAMEJUMP;
        end
        [shiftIdxS, shiftIdxE] = regexp(front, '(?<=_)[[0-9]]+(?=x)');
        Nshift = str2num(front(shiftIdxS(end):shiftIdxE(end)));
        if(rotation == 0)
            confName = sprintf('./%dx%d.conf.sk.png', Nshift, Nshift);
        else
            confName = sprintf('./%dx%d.conf.sk.rot%d.png', Nshift, Nshift, rotation);
        end

        groupCount = floor(imgCount / frameJump);
        if(Nshift ~= 1)
            if(groupMode == 1)
                groupCount = groupCount - Nshift * Nshift + 1;
            else
                groupCount = floor(groupCount / (Nshift * Nshift));
            end
        end
        clear centerviewList
        for startFrame = 0:frameJump - 1
            autoCenterFrame = startFrame;%自动找中心点用的那张图，默认取起始图片(zgx mode必须取扫描第一张），从0到N-1
            realignFolder = strcat(source{sourceIdx}(1:end-7), '__', num2str(startFrame));
            mkdir(realignFolder);
            realignFolderSub = strcat(realignFolder,'/realign');
            mkdir(realignFolderSub);
            realignName =  strcat(realignFolderSub, '/test');
            if(centerviewOnly == 1)
                centerView = strcat(realignFolder, '/centerView.only.tiff');%输出中心视角模式，若不需要，此处填一个'None'，需要则填tiff的名字， 若.tiff前有个.only，则只会输出中心视角（例如：outdir/a.only.tiff）
            else
                centerView = strcat(realignFolder, '/centerView.tiff');
            end

            centerviewList(startFrame+1, :) = strcat(centerView(1:end-4),'0.tif');%TODO: now we only support .0.tif

            if(autoCenterMode == 1)
                [centerX, centerY, centerProb] = AutoCenter(imread(name, autoCenterFrame + 1), Nnum);
            else
                centerProb = 999;
            end
            WriteCenterPt(realignFolder, centerX, centerY, centerProb, imread(name, autoCenterFrame + 1));
           
            % ./ReAlign [CaptureNshift] [OutputFile] [RawData]
            % [RawDataIdx] [CenterPixelX] [CenterPixelY]
            % [PickLensSizeX] (DataGroup = 1) (GroupMode = 0){0=individual, 1=continues} (ShiftConfigFile = Auto)
            % (FrameJump = 1) (Resize = 1) (ResizeShift = 13) (PickLensSizeY = X)
            command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d %s %s %d %f %f',...
            Nshift, realignName, name, startFrame, centerX, centerY, Nx,...
            groupCount, groupMode, confName, frameJump, resize, Nnum, Ny, realignMode, centerView, rotation, preResize, slightRotation);
            system(command);
        end
        if(frameJump > 1 && LaserMergeOption ~= 0)
            fprintf('working on channel merge...\n');
            channel = min([3, frameJump]);
            mergeView = uint16(zeros((Ny*2+1)*Nshift, (Nx*2+1)*Nshift, 3, groupCount));
            outPath = strcat(source{sourceIdx}(1:end-7), '__merge');
            mkdir(outPath);
            for frameIdx = 1:groupCount
                for channelIdx = 1:channel
                    if(exist(centerviewList(channelIdx, :), 'file'))
                        mergeView(:,:,channelIdx,frameIdx) = imread(centerviewList(channelIdx, :),  frameIdx);
                    end
                end
            end
            for channelIdx = 1:channel
                mergeView(:,:,channelIdx,:) = double(mergeView(:,:,channelIdx,:)) / double(max(max(max(max(mergeView(:,:,channelIdx,:)))))) * 255 * MergeContrastRatio;
            end
            mergeView = uint8(mergeView);
            for frameIdx = 1:groupCount
                if(frameIdx == 1)
                        imwrite(squeeze(mergeView(:,:,:,frameIdx)), strcat(outPath,'/centerViewMerge_channel.tif'));
                    else
                        imwrite(squeeze(mergeView(:,:,:,frameIdx)), strcat(outPath,'/centerViewMerge_channel.tif'), 'WriteMode', 'append');
                end
            end
        end
        fprintf('Realign %s, frame count : %d, cost %f sec\n', source{sourceIdx}, groupCount, toc(ticer));
    catch
       warning('we met some error!'); 
    end
end

function [] = WriteCenterPt(path, CX, CY, prob, img)
    txt = fopen(strcat(path,'/CenterPoint.txt'), 'w');
    fprintf(txt, 'X = %d\nY = %d\nprob = %.2f%%\n##Both (X,Y) in [0, N-1]', CX, CY, prob*100);
    fclose(txt);
    imwrite(img, strcat(path,'/UsedImg.tiff'));
end

function [Xcenter,Ycenter, prob] = AutoCenter(img, Nnum)

    SUB_NUM = 3;%n
    SUB_MAX_RANGE_RATIO = 5;%k
    MAX_AREA_TRUST_RATIO = 4;%m
    img = double(img);
    fullX = size(img, 2);
    fullY = size(img, 1);
    %整体，最亮处(kxk pix)(m倍可信度)，分块nxn
    ansList = zeros(1 + MAX_AREA_TRUST_RATIO + SUB_NUM*SUB_NUM, 3);
    fprintf('Full Image :');
    [ansList(1,1), ansList(1,2), ansList(1,3)] = AutoCenterSub(img, Nnum, 0, 0);
    fprintf('\n');
    maxImg = max(img(:));
    [maxPosY, maxPosX] = find(img == maxImg);
    maxPosX = maxPosX(round(size(maxPosX, 1) / 2));
    maxPosY = maxPosY(round(size(maxPosY, 1) / 2));
    
    rangeSX = round(maxPosX - fullX / SUB_MAX_RANGE_RATIO);
    if(rangeSX < 1) rangeSX = 1; end
    rangeEX = round(maxPosX + fullX / SUB_MAX_RANGE_RATIO);
    if(rangeEX > fullX) rangeEX = fullX; end
    
    rangeSY = round(maxPosY - fullY / SUB_MAX_RANGE_RATIO);
    if(rangeSY < 1) rangeSY = 1; end
    rangeEY = round(maxPosY + fullY / SUB_MAX_RANGE_RATIO);
    if(rangeEY > fullY) rangeEY = fullY; end
    fprintf('Lightest x%d:', MAX_AREA_TRUST_RATIO);
    [ansList(2,1), ansList(2,2), ansList(2,3)] = ...
        AutoCenterSub(img(rangeSY : rangeEY, rangeSX : rangeEX), Nnum, rangeSX - 1, rangeSY - 1);
    fprintf('\n');
    ansList(2:2+MAX_AREA_TRUST_RATIO-1, :) = repmat(ansList(2,:), [MAX_AREA_TRUST_RATIO, 1]);
    
    anchorPtX = zeros(SUB_NUM+1, 1);
    for i = 0:SUB_NUM
        anchorX = round(1 + i * fullX/SUB_NUM);
        if(anchorX > fullX) anchorX = fullX; end
        anchorPtX(i+1) = anchorX;
    end
    anchorPtY = zeros(SUB_NUM+1, 1);
    for i = 0:SUB_NUM
        anchorY = round(1 + i * fullY/SUB_NUM);
        if(anchorY > fullY) anchorY = fullY; end
        anchorPtY(i+1) = anchorY;
    end
    
    idx = 1 + MAX_AREA_TRUST_RATIO + 1;
    
    for i = 1:SUB_NUM%x
        for j = 1:SUB_NUM%y
            fprintf('Sub X=%d Y=%d:', i,j);
            [ansList(idx,1), ansList(idx,2), ansList(idx,3)] = ...
                AutoCenterSub(img(anchorPtY(j) : anchorPtY(j+1), anchorPtX(i) : anchorPtX(i+1)), Nnum, anchorPtX(i) - 1, anchorPtY(j) - 1);
            distance = sqrt((i - ((SUB_NUM + 1) /2))^2 + (j - ((SUB_NUM + 1) /2))^2);
            prob_loss = 1 - (distance / SUB_NUM);
            fprintf('prob loss ratio = %.2f\n',  prob_loss);
            ansList(idx,3) = ansList(idx,3) * prob_loss;
            idx = idx + 1;
        end
    end
    savedAnsList = ansList;
    ansList(:,3) = ansList(:,3) / sum(ansList(:,3));
    ansList(:,1) = ansList(:,1) .* ansList(:,3);
    ansList(:,2) = ansList(:,2) .* ansList(:,3);
    myAns = round([sum(ansList(:,1)), sum(ansList(:,2))]);
    prob = 0;
    probCount = 0;
    for i = 1:size(savedAnsList, 1)
        if(myAns == savedAnsList(i, 1:2))
            probCount = probCount + 1;
            prob = prob + savedAnsList(i, 3);
        end
    end
    if(probCount ~= 0)
        prob = prob / probCount;
    end
    Xcenter = myAns(1) + floor(size(img,2) / Nnum / 2) * Nnum;
    Ycenter = myAns(2) + floor(size(img,1) / Nnum / 2) * Nnum;
    fprintf('AutoCenter found x = %d, y = %d (range from 0~[size-1]), credibility = %.2f%%, check at ImageJ\n', Xcenter, Ycenter, prob*100);
end

%返回值，输入值（offset）均为0起点坐标
function [Xcenter, Ycenter, prob] = AutoCenterSub(img, Nnum, x_offset, y_offset)
    
    BEST_RATIO = 0.3;
    WORST_RATIO = 0.9;

    img = double(img);
    img = img ./ max(img(:));
    img = img.^2;
    kernal = fspecial('gaussian',[Nnum,Nnum],3);
    img = imfilter(img, kernal);
    locMatrix = zeros(Nnum, Nnum);
    for i = 1:Nnum
        for j = 1:Nnum
            picMat = img(i:Nnum:end, j:Nnum:end);
            %locMatrix(i,j) = mean(mean(picMat));
            avg = mean(mean(picMat));
            picMat(picMat < avg) = 0;
            hugePos = find(picMat ~= 0);
            locMatrix(i,j) = sum(picMat(:)) / size(hugePos, 1);
        end
    end
    
    sumX = sum(locMatrix);
    sumY = sum(locMatrix,2);
    sumX = sumX + circshift(sumX, 1, 2);
    sumY = sumY + circshift(sumY, 1, 1);
    darkX = 1;
    for i = 1:Nnum
        if(sumX(i) < sumX(darkX))
            darkX = i;
        end
    end
    darkY = 1;
    for i = 1:Nnum
        if(sumY(i) < sumY(darkY))
            darkY = i;
        end
    end
    Xcenter = mod(darkX + floor(Nnum / 2) - 1 + x_offset, Nnum);
    Ycenter = mod(darkY + floor(Nnum / 2) - 1 + y_offset, Nnum);
    prob = (WORST_RATIO - min(min(locMatrix)) / max(max(locMatrix))) / (WORST_RATIO - BEST_RATIO);
    if(prob > 1) 
        prob = 1; 
    elseif(prob < 0) 
        prob = 0; 
    end
    fprintf('AutoCenterSub x = %d, y = %d, prob = %f ', Xcenter, Ycenter, prob);
end



% function [CX, CY, prob] = AutoCenter(name, frame, Nnum)
% 
%     img = double(imread(name, frame + 1));
%     kernal = fspecial('gaussian',[Nnum,Nnum],3);
%     fed = imfilter(img, kernal);
%     locMatrix = zeros(Nnum, Nnum);
%     for i = 1:Nnum
%         for j = 1:Nnum
%             picMat = fed(i:Nnum:end, j:Nnum:end);
%             locMatrix(i,j) = mean(mean(picMat));
%         end
%     end
%     [yc,xc] = find(locMatrix == max(max(locMatrix)));
%     centerXX = xc + floor(size(img,2) / Nnum / 2) * Nnum - 1;
%     centerYY = yc + floor(size(img,1) / Nnum / 2) * Nnum - 1;
%     fprintf('AutoCenter found center point x = %d, y = %d (range from 0~[size-1]), check at ImageJ\n', centerXX, centerYY);
% %             if(mod(centerXX, Nnum) ~= mod(centerX, Nnum) || mod(centerYY, Nnum) ~= mod(centerY, Nnum))
% %                 warning('AutoCenter Point not match! Found(modded) = (%d, %d), input = (%d, %d), both at c++ range',...
% %                     mod(centerXX, Nnum), mod(centerYY, Nnum), mod(centerX, Nnum), mod(centerY, Nnum));
% %             end
%     locMatrixSort = sort(reshape(locMatrix,[1,Nnum * Nnum]));
%     yyt = 0; xxt = 0;
%     for i = 1:4
%         [yyyt, xxxt] = find(locMatrix == locMatrixSort(i));
%         yyt = yyyt / 4 + yyt;
%         xxt = xxxt / 4 + xxt;
%     end
%     %TODO:Here we are facing a problem, when 2 least pixel reach the edge!
%     prob = 0;
%     if(abs(abs(yyt - yc) - Nnum / 2) > 1e-3 || abs(abs(xxt - xc) - Nnum / 2) > 1e-3 || min(min(locMatrix)) / max(max(locMatrix)) > 0.9)
%         warning('呜呜呜，我找的中心点可能有问题……你能不能再检查一下用于找中心点的图呢QAQ');
%     else
%         ctemp = (0.9 - min(min(locMatrix)) / max(max(locMatrix))) / 0.2;
%         if(ctemp > 1) ctemp = 1; end
%         fprintf('AutoCenter calculation credibility = %.2f%%\n', ctemp * 100); %这个概率我瞎写的，大概就是最小值/最大值=0.9->0%（信噪比过低）;=0.7->100%（信噪比足够）
%         prob = ctemp;
%     end
%     CX = centerXX;
%     CY = centerYY;
%     %clear xc xxt xxxt yc yyt yyyt locMatrixSort img kernal fed ctemp
%     %t = imresize(locMatrix, [1300,1300], 'nearest');imshow(t,[]);
% 
% end
% 

