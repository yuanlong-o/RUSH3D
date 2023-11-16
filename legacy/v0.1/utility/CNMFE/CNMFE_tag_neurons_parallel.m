 function [obj, tags_] = CNMFE_tag_neurons_parallel(obj, min_pnr)
            if ~exist('min_pnr', 'var') || isempty(min_pnr)
                min_pnr = 3;
            end
            A_ = obj.A;
            %             C_ = obj.C_;
            S_ = obj.S;
            K = size(A_, 2);
            tags_ = zeros(K, 1, 'like', uint16(0));
            min_pixel = obj.options.min_pixel;
            
            % check the number of nonzero pixels
            nz_pixels = full(sum(A_>0, 1));
            tags_ = tags_ + uint16(nz_pixels'<min_pixel);
            
            % check the number of calcium transients after the first frame
            if obj.options.deconv_flag
                nz_spikes = full(sum(S_(:,2:end)>0, 2));
                tags_ = tags_ + uint16(nz_spikes<1)*2;
                
                tmp_std = std(obj.C_raw-obj.C, 0, 2);
                % check the noise level, it shouldn't be 0
                temp= tmp_std./GetSn(obj.C_raw);
                tags_ = tags_ + uint16(temp<0.1)*4;
                
                % check PNR of neural traces
                pnrs = max(obj.C, [], 2)./tmp_std;
                tags_ = tags_+uint16(pnrs<min_pnr)*8;
            end
            
            
            obj.tags = tags_;
        end
