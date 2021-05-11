function summary_image = preprocessing_module(realign_param, data_path, wdf_path)
%% this function conduct preprocessing of RUSH3D data.
%  last update: 5/11/2021. YZ


Nshift = realign_param.Nshift;
Nnum= realign_param.Nnum;

large_cycle = realign_param.large_cycle; %% 24 x 40 Images will be loaded;
small_cycle = realign_param.small_cycle; %% one stack contains 40 LF Images;
maxIter = realign_param.maxIter; %% Maximum times for Iterate rank normal one

%% Raw LF Data Path
%  path name illegal.
rawdata_folder1 = 'Z:/���ݲɼ�/';
rawdata_folder2 = 'gyd0309/';
rawdata_folder12 = strcat(rawdata_folder1, rawdata_folder2);
rawdatafilename_perfix = 'X1_NORMAL_3x3_25.0ms_HBar_Hardware_LaserCount1_210309202631';

%% WDF and STDLF Save Path

wdfdata_folder1 = 'D:/RichardW2021/RUSH3DResult/Preprocess/';
wdfdata_folder2 = rawdata_folder2;
wdfdata_folder12 = strcat(wdfdata_folder1, wdfdata_folder2);
wdfdata_folder3 = strcat(wdfdata_folder12, rawdatafilename_perfix,'/');
wdfdata_folder4 = strcat(wdfdata_folder3, 'WDF/');

if ~exist(wdfdata_folder12,'file')
    mkdir(wdfdata_folder12);
end

if ~exist(wdfdata_folder3,'file')
    mkdir(wdfdata_folder3);
end

if ~exist(wdfdata_folder4,'file')
    mkdir(wdfdata_folder4);
end

%% Preprocess + Save
name0 = [rawdatafilename_perfix,'.',num2str(0),'.tiff'];
name = [rawdata_folder12,name0];
tmp = double(imread(name,1));
tmp = imrotate(tmp,-0.029,'bicubic');


rawdata = zeros(size(tmp,1),size(tmp,2),(large_cycle+1)*small_cycle);
clear tmp;

m = 1;
for num = 0:1:large_cycle
    name0 = [rawdatafilename_perfix,'.',num2str(num),'.tiff'];
    name = [rawdata_folder12,name0];
    for i = 1:1:small_cycle
        tmp = double(imread(name,i));
        tmp = imrotate(tmp,-0.029,'bicubic');
        
        %tmp = tmp(cutx_start:cutx_end,cuty_start:cuty_end);
        
        rawdata(:,:,m) = tmp;
        if(mod(m,Nshift^2)==0)
            disp([num2str(m/(Nshift^2)), ' large has been loaded!']);
        end
        m = m+1;
    end
end

%% STD
meandata = zeros(size(tmp,1),size(tmp,2),Nshift^2);
for i = 1:Nshift^2
    tmp_stack = rawdata(:,:,i:Nshift^2:end);
    [bg_spatial, bg_temporal] = rank_1_NMF(reshape(tmp_stack,[size(tmp,1)*size(tmp,2),size(tmp_stack,3)]), maxIter);
    tmp_stack = tmp_stack - reshape(bg_spatial*bg_temporal,[size(tmp,1),size(tmp,2),size(tmp_stack,3)]);
    clear bg_spatial  bg_temporal;
    tmp_stack(tmp_stack<0) = 0;
    meandata(:,:,i) = sum(tmp_stack,3)./size(tmp_stack,3);
    if i == 1
        vardata = zeros(size(tmp,1),size(tmp,2),Nshift^2);
        stddata = zeros(size(tmp,1),size(tmp,2),Nshift^2);
    end
    vardata(:,:,i) =  sum((tmp_stack- repmat(meandata(:,:,i),[1,1,size(tmp_stack,3)])).^2,3);
    stddata(:,:,i) = sqrt(double(vardata(:,:,i))./(size(tmp_stack,3)-1));
    disp([num2str(i), ' std has been calculated!']);
end
imwriteTFSK(single(stddata),[wdfdata_folder3, 'rank_1_std.tif']);

%% Realign
tic
startFrame = 0; %��ʼ֡��(0~N-1)
frameJump = 1; %��֡����Զ�ɫ��ϣ����֡�������1~N)
stdfilename = [wdfdata_folder3,'rank_1_std.tif'];%����ĵ�һ��tiff�����֣����Զ������ң�����������ԭ�ļ�����
realignName = [wdfdata_folder4,'rank1_'];%���ͼ��ǰ׺����ע���Լ��½��ļ��У�����û����������ᱨ����
Nx = 120;%120ȡ��΢͸�������� Nx = (pick - 1) /2
Ny = 120;%ͬ�ϣ�Y���򣬴��������⣬����Nx==Ny
resize = 0;%�Ƿ�resize��13x13���
confName = ['./',num2str(Nshift),'x',num2str(Nshift),'.conf.sk.png'];%ɨ�������ļ�����Ҫ����Nshift���޸�
groupCount = 1;%���������
groupMode = 1;%0=���ģʽ(1-9,10-18,...);1=����ģʽ(1-9,2-10,...)
autoCenterMode = 0;%�Զ������ĵ�
autoCenterFrame = startFrame;%�Զ������ĵ��õ�����ͼ��Ĭ��ȡ��ʼͼƬ(zgx mode����ȡɨ���һ�ţ�����0��N-1
centerX=2079;%���ĵ�X���Զ�ģʽ��ʧЧ����0��N-1����c++����ϵ����0��ʼ����
centerY=2047;%���ĵ�Y���Զ�ģʽ��ʧЧ����0��N-1����c++����ϵ����0��ʼ����
realignMode = 'LZ';%realignMode = 'ZGX'; %LZ mode ���� ZGX mode
centerView = [wdfdata_folder3,'centerView_rank1.tiff'];%��������ӽ�ģʽ��������Ҫ���˴���һ��'None'����Ҫ����tiff�����֣� ��.tiffǰ�и�.only����ֻ����������ӽǣ����磺outdir/a.only.tiff��
rotation =  0;%��ת�ⳡͼ��˳ʱ�룬���0,90,180,270
preResize = 1;%Resize�ⳡͼ��Ԥ������Ĭ��ֵ1.0��
slightRotation = 0;%��΢��ת�ⳡͼ����ֵ�� �����΢��ת���������ת֮ǰ��������ֵ
%˳���ȸ�������ң������ĵ�����΢��ת����resize�ⳡͼ����resize���ĵ����꣬�ٲü�ͼ���ٽ��д���ת��Ȼ�����Realign��ȡ

if(autoCenterMode == 1)
    img = double(imread(stdfilename, autoCenterFrame + 1));
    [centerX, centerY, ~] = AutoCenter(img, Nnum);
end


command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d %s %s %d %f %f',...
    Nshift, realignName, stdfilename, startFrame, centerX, centerY, Nx,...
    groupCount, groupMode, confName, frameJump, resize, Nnum, Ny, realignMode, centerView, rotation, preResize, slightRotation);
system(command);
toc


% output the std file..
end

%% utility functions
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



function [Xcenter,Ycenter, prob] = AutoCenter(img, Nnum)

SUB_NUM = 3;%n
SUB_MAX_RANGE_RATIO = 5;%k
MAX_AREA_TRUST_RATIO = 4;%m
img = double(img);
fullX = size(img, 2);
fullY = size(img, 1);
%���壬������(kxk pix)(m�����Ŷ�)���ֿ�nxn
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

%����ֵ������ֵ��offset����Ϊ0�������