function  [raw_stack_after_rr, total_num] = rotate_resize_module_no_store(preprocess_param,input_rawdata_perfix,rr_rawdata_name)
%% Rotate_Resize_Module Rotate and resize the rawdata stack and save

% This program is used to rotate and resize raw data in one capture (one
% video), It loads all file (all stack) in sequence and save them in the
% same format after slightly rotate and resize. It return the result of one
% large stack at the end of the function mudule.

% last update: 6/5/2021. MW

% Input:
% preprocess_param        including total number of raw image(large_cycle*small_cycle) and rotate and resize parameter.
% input_rawdata_perfix
% rr_rawdata_name

% Output:
% total_num : total frame number
% output file             save in the outdir separately in several stacks
% disable the global stack since it can not be stored further.

% parameter parser
large_cycle = preprocess_param.large_cycle;
small_cycle = preprocess_param.small_cycle;
pre_rotate = preprocess_param.pre_rotate;
pre_resize = preprocess_param.pre_resize;

rawdata_name = strcat(input_rawdata_perfix,'.',num2str(0),'.tiff');
tmp = double(imread(rawdata_name,1));
tmp = imrotate(tmp,pre_rotate,'bicubic');
tmp = imresize(tmp,pre_resize,'bicubic');


m = 1;
for num = 0:1:large_cycle-1
    raw_stack_after_rr = zeros(size(tmp,1),size(tmp,2),small_cycle, 'single');
    if mod(num, 5) == 0
        fprintf('%d in %d processing \n', num, large_cycle-1)
    end
    for i = 1:1:small_cycle
        rawdata_name = strcat(input_rawdata_perfix,'.',num2str(num),'.tiff');
        tmp = double(imread(rawdata_name,i));
        % rotate
        tmp = imrotate(tmp,pre_rotate,'bicubic');
        % resize
        tmp = imresize(tmp,pre_resize,'bicubic');
        raw_stack_after_rr(:,:,i) = single(tmp);
        m = m+1;
    end
    imwriteTFSK(uint16(raw_stack_after_rr),[rr_rawdata_name,'.', num2str(num),'.tiff']);
    disp([num2str(m),' raw data has been rr and loaded']);
end
total_num = m - 1;
