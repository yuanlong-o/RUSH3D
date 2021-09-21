clc;clear;
addpath(genpath('utility'));
addpath('preprocess_module');
addpath('reconstruction_module');
addpath('registration');
debgrecon_param.mchannel = 1;
debgrecon_param.savegroup = 2;
preprocess_param.view_range = 5;
preprocess_param.Nnum = 15;
view_array = view_config(preprocess_param); 
[~, seq] = Spiral_circle(preprocess_param.Nnum,preprocess_param.view_range);
savegroup = debgrecon_param.savegroup;
gpuDevice;
%% Registration and std parameter
regstd_param.reg_gSig = 2;
regstd_param.reg_bin_width = 200;
regstd_param.rankdetrending = 0;
regstd_param.maxIter = 10;
maxframe = 1202;
%%
vid_debg_mul_video_savepath = 'D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test1\vid_mul_video'; 
outdir = 'D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test_reg_linear';
first_WDF = loadtiff('D:\RUSH3Dproject\RUSH3Dresult\0622\rasai148d_1_z287\test1\realign\realign_rasai148d_1_z287_3x3_45.0ms_Full_Hardware_LaserCount1_210622144634_No0.tif');
%% Registration and STD
disp('---------------------------Registration---------------------------');

reg_savepath = sprintf('%s\\reg_path', outdir);
if ~exist(reg_savepath,'file')
    mkdir(reg_savepath);
end

t_regstd_start = clock; 

std_WDF = zeros(size(first_WDF, 1), size(first_WDF, 2), length(view_array));
% for different channles
for ch_id = 1 : debgrecon_param.mchannel

    % measure shift
    tic;
    cv_ind = 113;
    centerview_video = [];
    for g_id = 1 : debgrecon_param.savegroup
        centerview_video = cat(3,centerview_video,double(loadtiff(sprintf('%s\\mul_video_view%d_ch%d_g%d.tiff',vid_debg_mul_video_savepath,cv_ind,ch_id,g_id))));
    end
    [d1,d2, ~] = size(centerview_video);
    [~, shifts, bound, option_r] = motion_correction(centerview_video, d1, d2, regstd_param, outdir);
    t_ms = toc;
    fprintf('measure shift process done and it takes %.2f secs\n',t_ms);
    %%
    % apply shift for each angle
    for v_ind = 1 : length(view_array)
        tic;
        v = find(v_ind == seq);
        view_video = [];
        for g_id = 1: debgrecon_param.savegroup
            view_video = cat(3,view_video,double(loadtiff(sprintf('%s\\mul_video_view%d_ch%d_g%d.tiff',vid_debg_mul_video_savepath, ...
                v,ch_id, g_id))));
        end
        t_s = clock;
        view_video = apply_shifts(view_video, shifts, option_r, bound/2, bound/2);
        t_se = etime(clock,t_s);
        disp(['shift: ',num2str(t_se), 's']);
        % save
        for g_id = 1: debgrecon_param.savegroup-1
            saveastiff(im2uint16(view_video(:,:,(g_id-1)*maxframe+1:g_id*maxframe)/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, v ,g_id));
        end
        saveastiff(im2uint16(view_video(:,:,(debgrecon_param.savegroup-1)*maxframe+1:end)/65535),sprintf('%s\\reg_view_%d_g_%d.tiff', reg_savepath, v ,savegroup));
        % calculate std
        t_std = clock;
        maxIter = regstd_param.maxIter;
        if regstd_param.rankdetrending == 1
            [curr_bg_spatial, curr_bg_temporal] = rank_1_NMF(reshape(view_video, [], size(view_video,3)), maxIter);
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)), ...
                curr_bg_spatial(:), curr_bg_temporal);
        else
            curr_std_image = compute_std_image(reshape(view_video, [], size(view_video,3)));
        end
        std_WDF(:, :, v) = reshape(curr_std_image, [size(view_video , 1), size(view_video , 2)]);
        t_stde = etime(clock,t_std);
        disp(['std: ',num2str(t_stde)]);
        t_onereg = toc;
        fprintf('%d in %d view has been registered and std and it takes %.2f secs\n', v, preprocess_param.Nnum^2, t_onereg);
    end
end
% save std
saveastiff(im2uint16(std_WDF / max(std_WDF(:))), sprintf('%s\\std_%s%d.tif', outdir, first_file_name,0));
t_regstd = etime(clock,t_regstd_start);
fprintf('resgistration process done and it takes %.2f secs\n',t_regstd);