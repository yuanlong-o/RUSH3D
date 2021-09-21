        function obj = CNMFE_remove_false_positives(obj, show_delete)
            if ~exist('show_delete', 'var')
                show_delete = false;
            end
            [obj, tags_] = CNMFE_tag_neurons_parallel(obj);  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
            ids = find(tags_);
            if isempty(ids)
                fprintf('all components are good \n');
            else
                if show_delete
                    obj = CNMFE_viewNeurons(obj, ids, obj.C_raw);
                else
                    obj = CNMFE_delete(obj, ids);
                end
            end
        end