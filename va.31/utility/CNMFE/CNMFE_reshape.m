        function Y = CNMFE_reshape(obj, Y, dim)
            % reshape the imaging data into diffrent dimensions
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            if dim==1
                Y=reshape(Y, d1*d2, []);  %each frame is a vector
            else
                Y = reshape(full(Y), d1, d2, []);    %each frame is an image
            end
        end