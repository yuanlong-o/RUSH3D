function [patch_info_array,wdf_h,wdf_w] = patch_preparation(first_WDF, patch_block, top_left, right_down)
% This function is used to segment the WDF into several patch at the beginning of the pipeline
% and save all the parameters in the patch info array 
first_WDF_cut = first_WDF(top_left(1):right_down(1),top_left(2):right_down(2));
wdf_h = size(first_WDF_cut,1);
wdf_w = size(first_WDF_cut,2);
patch_indx = round(linspace(top_left(1),right_down(1)+1,patch_block(1)+1));
patch_indy = round(linspace(top_left(2),right_down(2)+1,patch_block(2)+1));
patch_info_array = cell(prod(patch_block,'all'),1);
count = 1;
for i = 1 : patch_block(1)
    for j = 1 : patch_block(2)
        
        % loc
        patch_tl = [patch_indx(i),patch_indy(j)];
        patch_rd = [patch_indx(i+1)-1,patch_indy(j+1)-1];
        patch_size_h = patch_rd(1) - patch_tl(1) + 1;
        patch_size_w = patch_rd(2) - patch_tl(2) + 1;

        patch_loc = [patch_tl; patch_rd];
        patch_info.location = patch_loc;
        
        % size 
        patch_info.size = [patch_size_h, patch_size_w];
        
        % array
        patch_info_array{count} = patch_info;
        count = count + 1;
    end
end

end