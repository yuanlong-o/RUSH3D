function [ai, ind_success] = extract_a(HY, Y, ind_ctr, sz, spatial_constraints, min_pixel_val)
%% given a patch of raw & high-pass filtered calcium imaging data, extract
% spatial and temporal component of one neuron (ai, ci). if succeed, then
% return an indicator ind_succes with value 1; otherwise, 0.
%% inputs:
%       HY:     d X T matrix, filtered patch data. Note this is done in small patches!
%       Y:      d X T matrix, raw data
%       ind_ctr:        scalar, location of the center
%       sz:         2 X 1 vector, size of the patch
%       spatial_constraints: cell 
%% Author: Pengcheng Zhou, Carnegie Mellon University.

% modified by YZ. use relative intensity
%% parameters 
nr = sz(1);
nc = sz(2);
min_pixels = 5;

%% find pixels highly correlated with the center
% HY(HY<0) = 0;       % remove some negative signals from nearby neurons
tmp_corr = HY; 
% tmp_corr = tmp_corr / max(tmp_corr(:));
data = HY;
data(tmp_corr < (min_pixel_val * (max(tmp_corr(:)) - min(tmp_corr(:))) + min(tmp_corr(:)))) = 0;


%% estimate ci with the mean or rank-1 NMF
ci = mean(data, 1);

if norm(ci)==0  % avoid empty results 
    ai=[];
    ind_success=false;
    return;
end
%%%%%%%%% SIMPLER IS BETTER %%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate ai 
ai = data; 

% if spatial_constraints.circular
ai = circular_constraints(reshape(ai, nr, nc)); % assume neuron shapes are spatially convex
% end

if spatial_constraints.connected
    ai = connectivity_constraint(reshape(ai, nr, nc));
end
ai = ai(:); 

% %% threshold the spatial shape and remove outliers 
% % remove outliers 
% temp =  full(ai>quantile(ai(:), 0.5)); 
% l = bwlabel(reshape(temp, nr, nc), 4); 
% temp(l~=l(ind_ctr)) = false; 
% ai(~temp(:)) = 0; 
if sum(ai(:)>0) < min_pixels %the ROI is too small
    ind_success=false;
    return;
end

% refine ci given ai 
% ind_nonzero = (ai>0);
% ai_mask = mean(ai(ind_nonzero))*ind_nonzero;
% ci0 = (ai-ai_mask)'*ai\((ai-ai_mask)'*Y);
% plot(ci, 'r'); 
% pause; 

% we use two methods for estimating the noise level 

% ind_neg = (ci<-4*sn); 
% ci(ind_neg) = rand(sum(ind_neg), 1)*sn; 

% normalize the result
% ci = ci / sn;
% ai = ai * sn;
% % return results
if norm(ai)==0
    ind_success= false;
else
    ind_success=true;
end
