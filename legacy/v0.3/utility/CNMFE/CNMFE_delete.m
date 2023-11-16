 function obj = CNMFE_delete(obj, ind)
            % write the deletion into the log file
            if ~exist('ind', 'var') || isempty(ind)
                return;
            end
            obj.A(:, ind) = [];
            obj.C(ind, :) = [];
            if ~isempty(obj.S)
                try obj.S(ind, :) = []; catch; end
            end
            if ~isempty(obj.C_raw)
                try obj.C_raw(ind, :) = []; catch;  end
            end
            if isfield(obj.P, 'kernel_pars')&&(  ~isempty(obj.P.kernel_pars))
                try obj.P.kernel_pars(ind, :) = []; catch; end
            end
            try  obj.ids(ind) = [];   catch;   end
            try obj.tags(ind) =[]; catch; end
            
            % save the log
        end