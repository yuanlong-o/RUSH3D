function Ybg = CNMFE_reconstruct_background(obj, Y_mat, frame_range)
            %%reconstruct background using the saved data
            % input:
            %   frame_range:  [frame_start, frame_end], the range of frames to be loaded
            %% Author: Pengcheng Zhou, Columbia University, 2017
            %% email: zhoupc1988@gmail.com
            
            %% process parameters
            
            try
                % map data
                d1 = size(Y_mat, 1);
                d2 = size(Y_mat, 2);
                T = size(Y_mat, 3);
                obj.options.d1 = d1;
                obj.options.d2 = d2;

                patch_pos = [[1; d1], [1; d2]];    % patch position
                block_pos = [[1; d1], [1; d2]];    % patch position
            catch
                error('No data file selected');
            end
            
            if ~exist('frame_range', 'var')||isempty(frame_range)
                frame_range = obj.frame_range;
            end
            % frames to be loaded for initialization
            T = diff(frame_range) + 1;
            
            bg_model = obj.options.background_model;
            bg_ssub = obj.options.bg_ssub;
            % reconstruct the constant baseline
            if strcmpi(bg_model, 'ring')
                b0_ = reconstruct_b0(Y_mat, obj);
                b0_new_ = CNMFE_reshape(obj, obj.b0_new, 2);
            end
            
            %% start updating the background
            Ybg = zeros(d1, d2, T);

            if strcmpi(bg_model, 'ring')
                W_ring = obj.W;
                %                     b0_ring = obj.b0;
                % load data
                Ypatch = Y_mat(:, :, frame_range(1) : frame_range(2));
                [nr_block, nc_block, ~] = size(Ypatch);
                Ypatch = reshape(Ypatch, [], T);
                tmp_block = block_pos;
                tmp_patch = patch_pos;
                b0_ring = b0_(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4));
                b0_ring = reshape(b0_ring, [], 1);

                b0_patch = reshape(b0_new_(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)), [], 1);

                % find the neurons that are within the block
                mask = zeros(d1, d2);
                mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
                ind = (reshape(mask(:), 1, [])* obj.A_prev>0);

                A_patch = obj.A_prev(logical(mask), ind);
                C_patch = obj.C_prev(ind, frame_range(1):frame_range(2));

                % reconstruct background
                %                     Cmean = mean(C_patch , 2);
                Ypatch = bsxfun(@minus, double(Ypatch), b0_ring);
                %                     b0_ring = b0_(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4));
                %                     b0_ring = reshape(b0_ring, [], 1);
                %
                if bg_ssub==1
                    Bf = W_ring*(double(Ypatch) - A_patch*C_patch);
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, Bf, b0_patch), diff(tmp_patch(1:2))+1, [], T);
                else
                    [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                    temp = reshape(double(Ypatch)-A_patch*C_patch, nr_block, nc_block, []);
                    temp = imresize(temp, 1./bg_ssub, 'nearest');
                    Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                    Bf = imresize(Bf, [nr_block, nc_block], 'nearest');
                    Bf = Bf((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1, :);
                    Bf = reshape(Bf, [], T);
                    Ybg(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4),:) = reshape(bsxfun(@plus, Bf, b0_patch), diff(tmp_patch(1:2))+1, [], T);
                end
            else
                error('only ring mode supported!')
            end

        end
   function b0_ = reconstruct_b0(Y_crop, obj)
    try
        % dimension of data
        d1 = size(Y_crop, 1);
        d2 = size(Y_crop, 2);
        T = size(Y_crop, 3);
        obj.options.d1 = d1;
        obj.options.d2 = d2;

        % parameters for patching information
        patch_pos =[[1; d1], [1; d2]];
        % number of patches

    catch
        error('No data file selected');
        b0_= [];
        return;
    end

    b0_ = zeros(d1, d2);
    b0_patch = obj.b0;
    tmp_patch = patch_pos;
    r0 = tmp_patch(1);
    r1 = tmp_patch(2);
    c0 = tmp_patch(3);
    c1 = tmp_patch(4);
	b0_(r0:r1, c0:c1) = reshape(b0_patch, r1-r0+1, c1-c0+1);

end
