function [obj, ind_del] = CNMFE_viewNeurons(obj, ind, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A~=0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = CNMFE_estCenter(obj);      %neuron's center
gSiz = obj.options.gSiz;        % maximum size of a neuron

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end


%% start viewing neurons
figure('position', [100, 100, 1024, 512]);
m=1;
while and(m>=1, m<=length(ind))
    %% full-frame view
    subplot(221); cla;
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    axis equal; axis off;
    if ind_del(m)
        title(sprintf('Neuron %d', ind(m)), 'color', 'r');
    else
        title(sprintf('Neuron %d', ind(m)));
    end
    %% zoomed-in view
    subplot(222); cla;
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end
    
    %% temporal components
    subplot(2,2,3:4);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m))), 'linewidth', 2); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    xlim([t(1), t(end)]);
    xlabel(str_xlabel);
    
    %% save images
    if save_img
        drawnow();
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d, keep(k, default)/delete(d)/split(s)/trim(t)\n\t/trim cancel(tc)/delete all(da)/backward(b)/end(e):    ', ind(m));
        
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'da')
            ind_del(m:end) = true;
            break;
        elseif strcmpi(temp, 'k')
            ind_del(m) = false;
            m= m+1;
        elseif strcmpi(temp, 's')
            try
                subplot(222);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                tmpA = obj.A(:, ind(m));
                obj.A(:, end+1) = tmpA.*tmp_ind(:);
                obj.C(end+1, :) = obj.C(ind(m), :);
                obj.A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                obj.S(end+1, :) = obj.S(ind(m), :);
                obj.C_raw(end+1, :) = obj.C_raw(ind(m), :);
                obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(m), :);
                k_ids = obj.P.k_ids;
                obj.ids(end+1) = k_ids+1;   % assign an neuron id
                obj.tags(end+1) = obj.tags(ind(m));
                obj.P.k_ids = k_ids+1;
            catch
                fprintf('the neuron was not split\n');
            end
        elseif strcmpi(temp, 't')
            try
                subplot(222);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                Amask(:, ind(m)) = tmp_ind(:);
                ind_trim(m) = true;
            catch
                fprintf('the neuron was not trimmed\n');
            end
        elseif strcmpi(temp, 'tc')
            Amask(:, ind(m)) = (obj.A(:, ind(m)) > 0);
            ind_trim(m) = false;
        elseif strcmpi(temp, 'e')
            break;
        elseif ~isnan(str2double(temp))
            m = m + floor(str2double(temp));
            m = max(m, 1);
            m = min(m, length(ind));
            fprintf('jump to neuron %d / %d\n', m, length(ind));
        else
            m = m+1;
        end
    end
end
if save_img
    cd(cur_cd);
else
    if ~isempty(ind(ind_trim))
        obj.A(:, ind(ind_trim)) = obj.A(:,ind(ind_trim)).*Amask(:, ind(ind_trim));
    end
    
    if ~isempty(ind(ind_del))
        obj = CNMFE_delete(obj, ind(ind_del));
    end
    %     obj.Coor = obj.get_contours(0.9);
      
    return;
end
end

