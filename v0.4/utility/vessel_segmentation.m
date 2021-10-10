function response_stack_segm = vessel_segmentation(img,outdir)

symmfilter = struct();
symmfilter.sigma     = 8; % variance of DoG %8
symmfilter.len       = 100; % rho, querying size %100
symmfilter.sigma0    = 1; % blurress kernel %1
symmfilter.alpha     = 0.15; % distance decreasing % 0.2
asymmfilter = false;

img = max(img(:)) - img;
response_stack = BCOSFIRE_lfm(img, symmfilter, asymmfilter); % only 2d here
% figure, imshow(response_stack)

% calculate the mask
se = strel('disk',1);
response_stack_segm= imdilate(response_stack, se);
response_stack_segm = response_stack_segm > 5e-2; %5
% response_stack_segm = imgaussfilt(response_stack_segm, 10);
% figure, imshow(response_stack_segm)
saveastiff(im2uint16(response_stack_segm), sprintf('%s\\blood_vessel_mask.tiff', outdir));