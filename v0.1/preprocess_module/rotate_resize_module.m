function [raw_stack_after_rr] = rotate_resize_module(preprocess_param,input_rawdata_perfix,rr_rawdata_name)
%% Rotate_Resize_Module Rotate and resize the rawdata stack and save

% This program is used to rotate and resize raw data in one capture (one
% video), It loads all file (all stack) in sequence and save them in the
% same format after slightly rotate and resize. It return the result of one
% large stack at the end of the function mudule.

% last update: 05/15/2021. MW

% Input:
% preprocess_param        including total number of raw image(large_cycle*small_cycle) and rotate and resize parameter.
% input_rawdata_perfix
% rr_rawdata_name

% Output:
% raw_stack_after_rr
% output file             save in the outdir separately in several stacks

large_cycle = preprocess_param.large_cycle;
small_cycle = preprocess_param.small_cycle;
pre_rotate = preprocess_param.pre_rotate;
pre_resize = preprocess_param.pre_resize;

rawdata_name = strcat(input_rawdata_perfix,'.',num2str(0),'.tiff');
tmp = double(imread(rawdata_name,1));
tmp = imrotate(tmp,pre_rotate,'bicubic');
tmp = imresize(tmp,pre_resize,'bicubic');
raw_stack_after_rr = zeros(size(tmp,1),size(tmp,2),large_cycle*small_cycle);

m = 1;
for num = 0:1:large_cycle-1
    for i = 1:1:small_cycle
        rawdata_name = strcat(input_rawdata_perfix,'.',num2str(num),'.tiff');
        tmp = double(imread(rawdata_name,i));
        tmp = imrotate(tmp,pre_rotate,'bicubic');
        tmp = imresize(tmp,pre_resize,'bicubic');
        raw_stack_after_rr(:,:,m) = tmp;
        m = m+1;
    end
    imwriteTFSK(uint16(raw_stack_after_rr(:,:,m - small_cycle:m-1)),[rr_rawdata_name,'.', num2str(num),'.tiff']);
    disp([num2str(m),' raw data has been rr and loaded']);
end