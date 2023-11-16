function [img, col0, AA] = CNMFE_overlapA(obj, ind, ratio)
            %merge all neurons' spatial components into one singal image
            if nargin<2 || isempty(ind)
                AA = obj.A;
            else
                AA = obj.A(:, ind);
            end
            if nargin<3
                ratio = 0.3;
            end
            
            %             v_max = max(max(AA,1));
            v_max = 1;
            AA = bsxfun(@times, AA, v_max./max(AA,[],1));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            [d, K] = size(AA);
            
            if K==2
                col = [4, 2];
            elseif K==3
                col = [4, 2, 1];
            else
                col = randi(6, 1, K);
            end
            col0 = col;
            img = zeros(d, 3);
            for m=3:-1:1
                img(:, m) = sum(bsxfun(@times, AA, mod(col, 2)), 2);
                col = floor(col/2);
            end
            img = obj.reshape(img, 2);
            img = img/max(img(:))*(2^16);
            img = uint16(img);
        end