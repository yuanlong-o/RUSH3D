function sn = CNMFE_estimate_noise(Y_crop, frame_range, method)

            d1 = size(Y_crop, 1);
            d2 = size(Y_crop, 2);
            T = size(Y_crop, 3);
            if ~exist('frame_range', 'var') || isempty(frame_range)
                frame_range = [1, min(T, 3000)];
            end
            T = diff(frame_range)+1;
            if ~exist('method', 'var') || isempty(method)
                method = 'psd';
            end
            if strcmpi(method, 'hist')
                % fit the histogram of with a parametric gaussian shape.
                % it works badly when the
                % most frames are involved in calcium transients.
                foo = @(Y) GetSn_hist(Y, false);
            elseif strcmpi(method, 'std')
                % simly use the standard deviation. works badly when
                % neurons have large calcium transients
                foo = @(Y) std(Y, 0, 2);
            else
                % default options is using the power spectrum method.
                % works badly when noises have temporal correlations
                foo = @CNMFE_GetSn;
            end
          
            
            Ypatch = Y_crop;
            Ypatch = double(reshape(Ypatch, [], T));
            tmp_sn = reshape(foo(Ypatch), d1, d2);

%             fprintf('Time cost for estimating the nosie levels:  %.3f \n\n', toc);
            sn = tmp_sn;
        end