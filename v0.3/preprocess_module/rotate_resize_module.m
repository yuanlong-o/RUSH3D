function  [raw_stack_after_rr, total_num] = rotate_resize_module(preprocess_param,input_rawdata_perfix,rr_rawdata_name,bg_rawdata_name)
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
% total_num : total frame number
% output file             save in the outdir separately in several stacks

% parameter parser
Nshift = preprocess_param.Nshift;
large_cycle = preprocess_param.large_cycle;
small_cycle = preprocess_param.small_cycle;
pre_rotate = preprocess_param.pre_rotate;
pre_resize = preprocess_param.pre_resize;
subnoise_flag = preprocess_param.subnoise_flag;

rawdata_name = strcat(input_rawdata_perfix,'.',num2str(0),'.tiff');
tmp = double(imread(rawdata_name,1));
tmp = imrotate(tmp,pre_rotate,'bicubic');
tmp = imresize(tmp,pre_resize,'bicubic');
raw_stack_after_rr = zeros(size(tmp,1),size(tmp,2),large_cycle*small_cycle, 'single');

if subnoise_flag == 1
    bg = double(imread(strcat(bg_rawdata_name,'.tiff'),1));
    bg_stack = zeros(size(bg,1),size(bg,2),Nshift^2);
    for i = 1: Nshift^2
        bg_stack(:,:,i) = double(imread(strcat(bg_rawdata_name,'.tiff'),i));
    end
end

m = 1;
if subnoise_flag == 1
    % sub ground noise
    for num = 0:1:large_cycle-1
        if mod(num, 5) == 0
            fprintf('%d in %d processing \n', num, large_cycle-1)
        end
        for i = 1:1:small_cycle
            rawdata_name = strcat(input_rawdata_perfix,'.',num2str(num),'.tiff');
            tmp = double(imread(rawdata_name,i));
            tmp = max(tmp - bg_stack(:,:,mod((num*large_cycle+i-1),Nshift^2)+1),0);
            tmp = imrotate(tmp,pre_rotate,'bicubic');
            tmp = imresize(tmp,pre_resize,'bicubic');
            raw_stack_after_rr(:,:,m) = single(tmp);
            m = m+1;
        end
        imwriteTFSK(uint16(raw_stack_after_rr(:,:,m - small_cycle:m-1)),[rr_rawdata_name,'.', num2str(num),'.tiff']);
        disp([num2str(m-1),' raw data has been rr and loaded']);
    end
    
else
    
    % do not sub ground noise
    for num = 0:1:large_cycle-1
        if mod(num, 5) == 0
            fprintf('%d in %d processing \n', num, large_cycle-1)
        end
        for i = 1:1:small_cycle
            rawdata_name = strcat(input_rawdata_perfix,'.',num2str(num),'.tiff');
            tmp = double(imread(rawdata_name,i));

            tmp = imrotate(tmp,pre_rotate,'bicubic');
            tmp = imresize(tmp,pre_resize,'bicubic');
            raw_stack_after_rr(:,:,m) = single(tmp);
            m = m+1;
        end
        imwriteTFSK(uint16(raw_stack_after_rr(:,:,m - small_cycle:m-1)),[rr_rawdata_name,'.', num2str(num),'.tiff']);
        disp([num2str(m-1),' raw data has been rr and loaded']);
    end
    
end
total_num = m - 1;
