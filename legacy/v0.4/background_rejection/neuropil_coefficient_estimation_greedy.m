function [coefs, auto_trace, ca_out, sp_out] = neuropil_coefficient_estimation_greedy(center_trace, ...
    neuropil_trace, lam, gama, lambda_l0)
%% this approach tries the naive coefficient estimation via greedy method.
% last update: 4/24.

% lambda controls the OASIS approach
% lambda_l0 control the lo regularizer (recommand 1e-2 when lam is 0.5)

%% pre justify
max_center_trace = max(center_trace);
center_trace = center_trace / max_center_trace;
center_trace = center_trace(:);
neuropil_trace = neuropil_trace / max(neuropil_trace)';
neuropil_trace  = neuropil_trace(:);


% OASIS
maxNeurop = 1.5;
decimate = 1;
% gama = 0.9;
% lam = 0.5;
maxIter_oasis = 100;

% lambda_l0 = 1e-2;
%% 
bin_num = 50;
coefs_array = linspace(0, maxNeurop, bin_num);
% close all
for k = 1:bin_num
    Fsub = center_trace - coefs_array(k) * neuropil_trace; % do subtraction
    
    % OASIS
    [ca, sp, Fbase, g, active_set] = foopsi_oasisAR1_ext(Fsub, gama, lam, false, ...
        true, decimate, maxIter_oasis);
    % residual
    res_deconv = Fsub - ca;
    data_fidelity = norm(res_deconv);
    %  spike number 
    spike_n = sum(sp>0);
    regularizer = spike_n;
    %  spike amplitude
    
    % record loss
    
    loss(k) = lambda_l0 *regularizer +  data_fidelity;
end
% find min
[~, min_idx] = min(loss);
coefs = coefs_array(min_idx);
auto_trace = center_trace - coefs * neuropil_trace;

% figure(105), plot(coefs_array, loss, 'linewidth', 1.5), xlabel('coefficient'), ylabel('loss')
% set(gca, 'Fontsize', 18)
% 
% figure(102), plot(center_trace - coefs * neuropil_trace, 'm', 'linewidth', 1.5)
[ca_out, sp_out, ~, ~, ~] = foopsi_oasisAR1_ext(auto_trace, gama, lam, false, ...
    true, decimate, maxIter_oasis);

auto_trace = max_center_trace * auto_trace;
ca_out = max_center_trace * ca_out;
% hold on, plot(ca, 'b', 'linewidth', 1.5)
% hold on, plot(sp, 'k', 'linewidth', 1.5')
% 
% % manual
% manual_adj = 1;
% figure(101), plot(center_trace + 2, 'r', 'linewidth', 1.5);
% hold on, plot(neuropil_trace + 1, 'g', 'linewidth', 1.5);
% hold on, plot(center_trace - manual_adj * neuropil_trace, 'b', 'linewidth', 1.5)
% hold on, plot(center_trace - coefs * neuropil_trace, '--m', 'linewidth', 1.5)
% 