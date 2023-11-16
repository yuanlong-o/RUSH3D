function [Mpr, shifts, bound, option_nr] = motion_correction(Yf, d1, d2, regstd_param, outdir)
%% this file use modules from normcorre to correct both burst motion and slow
%  motion of the input movie.


%  last update: 6/17/2020. YZ

%% pre processing
gSig = regstd_param.reg_gSig; 
gSiz = 3*gSig;
bin_width = regstd_param.reg_bin_width;

psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
%Y = imfilter(Yf,psf,'same');
%bound = 2*ceil(gSiz/2);
Y = imfilter(Yf,psf,'symmetric');
bound = 0;

%% parameters
% parsed parameter
option_r.d1 = d1-bound;
option_r.d2 = d2 - bound;
option_r.d3 = 1;
option_r.bin_width = bin_width; % 


% default parameter
option_r.grid_size = [option_r.d1, option_r.d2, 1];   % grid size
option_r.overlap_pre = [20, 20, 1]; % size of overlapped region before upsampling
option_r.overlap_post = [20, 20, 1]; % size of overlapped region after upsampling
option_r.min_patch_size = [10,10,16];%[32, 32, 16]; % minimize patch size
option_r.min_diff = [16, 16, 5];
option_r.us_fac = 50; % upsampling factor
option_r.mot_uf = [1, 1, 1];  % degree of patches upsampling (default: [4,4,1])
option_r.max_dev = [3, 3, 1];
option_r.max_shift = [20, 20, 20];
option_r.init_batch = 100;
option_r.buffer_width = 50;
option_r.iter = 3;

option_r.init_batch = 100; % length of init batch

option_r.window_length = 0.5;
option_r.nFrames = 50; % frames to avg

option_r.mem_batch_size = 1000;
option_r.add_value = 0; % flat, add dc value to data

option_r.upd_template = true; %  flag for online template updating
option_r.phase_flag = false;
option_r.print_msg = false; % flat for printing
option_r.boundary = 'copy';
option_r.shifts_method = 'fft';
option_r.correct_bidir = false; 
option_r.use_windowing = false;
option_r.output_type = 'mat';
option_r.method = {'median', 'mean'}; % average the template, median or mean
option_r.col_shift = [];
option_r.memmap = false;
option_r.mem_filename = [];

buf = size(Yf);
T = buf(end);
%% rigid shift
% register, this is quite quick. Output would be a scalar
tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),...
    option_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
tic; Mr = apply_shifts(Yf,shifts1,option_r,bound/2,bound/2); toc % apply shifts to full dataset

% metrics, on filtered image
[cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),option_r.max_shift);
[cM1,mM1,vM1] = motion_metrics(M1,option_r.max_shift);
[cYf,mYf,vYf] = motion_metrics(Yf,option_r.max_shift); % Yf maybe stands for Y full?
[cM1f,mM1f,vM1f] = motion_metrics(Mr,option_r.max_shift);

% figure;
% shifts_r = squeeze(cat(3,shifts1(:).shifts));
% subplot(311); plot(shifts_r);
%     title('Rigid shifts','fontsize',14,'fontweight','bold');
%     legend('y-shifts','x-shifts');
% subplot(312); plot(1:T,cY,1:T,cM1);
%     title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
%     legend('raw','rigid');
% subplot(313); plot(1:T,cYf,1:T,cM1f);
%     title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
%     legend('raw','rigid');
% print(fullfile(outdir, [datestr(now, 'YYmmddTHHMM') '_rigid_registration.pdf']), '-dpdf', '-r300');
%% non-rigid shift
option_nr = option_r;
option_nr.bin_width = round(bin_width / 2); % better resolution
option_nr.grid_size = [min(floor(d1/2),256),min(floor(d1/2),256), 1];%256
option_nr.mot_uf = [4, 4, 1];
option_nr.correct_bidir = false;
option_nr.overlap_pre = [32, 32, 1];%32
option_nr.overlap_post =[32, 32, 1];%32
option_nr.iter = 3;



tic; [M2,shifts,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:), ...
    option_nr, template1); toc % register filtered data, supply rigid template
tic; Mpr = apply_shifts(Yf,shifts, option_nr, bound/2, bound/2); toc % apply the shifts to the removed percentile


% metrics
[cM2,mM2,vM2] = motion_metrics(M2,option_nr.max_shift);

%% plo

% shift
shifts_r = squeeze(cat(3,shifts1(:).shifts)); % rigid shift
shifts_nr = cat(ndims(shifts(1).shifts)+1,shifts(:).shifts);% non-rigid shift
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,2,:))';
shifts_y = squeeze(shifts_nr(:,1,:))';


% figure;
% ax1 = subplot(311);
% plot(1:T, cY,1:T, cM1, 1:T,cM2); 
% legend('raw data','rigid','non-rigid'); title('correlation coefficients for filtered data','fontsize',14,'fontweight','bold') 
% set(gca,'Xtick',[],'XLim',[0,T-3])
% 
% ax2 = subplot(312); 
% plot(shifts_x); 
% hold on; plot(shifts_r(:,2),'--k','linewidth',2); 
% title('displacements along x','fontsize',14,'fontweight','bold')
% set(gca,'Xtick',[])
% 
% ax3 = subplot(313); 
% plot(shifts_y); 
% hold on; plot(shifts_r(:,1),'--k','linewidth',2); 
% 
% title('displacements along y','fontsize',14,'fontweight','bold')
% xlabel('timestep','fontsize',14,'fontweight','bold')
% linkaxes([ax1,ax2,ax3],'x')
% 
% print(fullfile(outdir, [datestr(now, 'YYmmddTHHMM') '_non_rigid_registration.pdf']), '-dpdf', '-r300');
end