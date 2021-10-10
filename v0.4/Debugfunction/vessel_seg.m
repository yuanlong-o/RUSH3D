outdir = 'D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test2';
std_WDF = loadtiff(sprintf('%s\\realign\\realign_rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_210622144634_No0.tif',outdir));
size_v = size(std_WDF,3);
img = double(std_WDF(:,:,ceil(size_v/2)));
% img = imresize(img,1/3);
img = double(img/max(img(:)));
figure,imshow(img,[]);
hold on;
%% blood vessel mask
% run blood vessel extraction
symmfilter = struct();
symmfilter.sigma     = 8; % variance of DoG %8
symmfilter.len       = 100; % rho, querying size %100
symmfilter.sigma0    = 1; % blurress kernel %1
symmfilter.alpha     = 0.25; % distance decreasing % 0.2
asymmfilter = false;

% down sample the image
% image_d = imresize(image, 0.2);

% Apple BCOSFIRE filter
% img = movie(:, :, 1);
img = max(img(:)) - img;
response_stack = BCOSFIRE_lfm(img, symmfilter, asymmfilter); % only 2d here
figure, imshow(response_stack)

% calculate the mask
se = strel('disk',1);
response_stack_segm= imdilate(response_stack, se);
response_stack_segm = response_stack_segm > 5e-2; %5
%%
% response_stack_segm = imgaussfilt(response_stack_segm, 10);
figure, imshow(response_stack_segm)
figure, imshow((1-response_stack_segm).*img,[]);
saveastiff(im2uint16(response_stack_segm), sprintf('%s\\blood_vessel_mask.tiff', outdir));