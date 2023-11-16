        function [center] = CNMFE_estCenter(obj)
            center = com(obj.A, obj.options.d1, obj.options.d2);
        end
        