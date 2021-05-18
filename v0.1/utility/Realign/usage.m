
tic
startFrame = 0; %起始帧号(0~N-1)
frameJump = 1; %跳帧：针对多色和希望跳帧的情况（1~N)
name = 'mouse_liver_8_9hz_3x3_100.0ms_Full_Hardware_LaserCount1_190903213802.0.tiff';%输入的第一组tiff的名字（会自动往下找，请勿重命名原文件名）
realignName = 'realign/test';%输出图像前缀（请注意自己新建文件夹，否则没有输出，不会报错）
mkdir('realign');%新建输出的文件夹
Nshift = 5;%拍摄的扫描
Nx = 75;%取的微透镜阵列数 Nx = (pick - 1) /2
Ny = 75;%同上，Y方向，代码有问题，建议Nx==Ny
Nnum = ;%微透镜阵列后的pixel数，一般不需要动
resize = 0;%是否resize到13x13输出
confName = './3x3.conf.sk.png';%扫描配置文件，需要根据Nshift做修改
groupCount = 220;%计算多少组
groupMode = 0;%0=间隔模式(1-9,10-18,...);1=连续模式(1-9,2-10,...)
autoCenterMode = 0;%自动找中心点
autoCenterFrame = startFrame;%自动找中心点用的那张图，默认取起始图片(zgx mode必须取扫描第一张），从0到N-1
centerX=1131;%中心点X，自动模式下失效，从0到N-1（在c++坐标系，从0开始！）
centerY=1087;%中心点Y，自动模式下失效，从0到N-1（在c++坐标系，从0开始！）
realignMode = 'LZ';%realignMode = 'ZGX'; %LZ mode 或者 ZGX mode
centerView = 'realign/centerView.only.tiff';%输出中心视角模式，若不需要，此处填一个'None'，需要则填tiff的名字， 若.tiff前有个.only，则只会输出中心视角（例如：outdir/a.only.tiff）
rotation =  0;%旋转光场图，顺时针，可填：0,90,180,270
preResize = 1;%Resize光场图（预处理，默认值1.0）
slightRotation = 0;%轻微旋转光场图（插值） 这个轻微旋转在上面的旋转之前做，带插值
%顺序：先根据所填（找）的中心点做轻微旋转，再resize光场图，再resize中心点坐标，再裁剪图像，再进行大旋转，然后进入Realign抽取

if(autoCenterMode == 1)
    img = double(imread(name, autoCenterFrame + 1));
    [centerX, centerY, ~] = AutoCenter(img, Nnum);
end


command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d %s %s %d %f %f',...
Nshift, realignName, name, startFrame, centerX, centerY, Nx,...
groupCount, groupMode, confName, frameJump, resize, Nnum, Ny, realignMode, centerView, rotation, preResize, slightRotation);
system(command);
toc


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