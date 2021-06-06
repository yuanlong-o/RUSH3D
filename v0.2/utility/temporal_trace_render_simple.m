function temporal_trace_render_simple(C)

ts = zscore(C , 0, 2);
% ts = zscore(global_T_deconv_filtered( 10500 : 10550, 1 : end) , 0, 2);
for i = 1 : size(ts, 1)
    ts(i, :) = smooth(ts(i, :).', 3);
end
y_shift = 2;
clip = false;
% rng(10021)
sel = 1:size(ts,1);

nixs = 1:size(ts,1);
sel_nixs = nixs(sel);
hold on
for n_ix = 1:numel(sel_nixs)
    ax = gca();
    ax.ColorOrderIndex = 1;
    loop_ts = ts(sel_nixs(n_ix),:);
    if clip
        loop_ts(loop_ts > 3*y_shift) = y_shift;
        loop_ts(loop_ts < -3*y_shift) = -y_shift;
    end
    t = (0:size(ts,2)-1);
    if max(loop_ts) - mean(loop_ts) > 2.5 * y_shift
        
        plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'color', 'k');
        plot1.Color(4) = 0.5;
    else
        % with alpha
        plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'k');
        plot1.Color(4) = 0.3;
    end
    
%     text(30, y_shift*(n_ix-1), num2str(sel_nixs(n_ix)));
end
xlabel('Frame');
xlim([min(t) max(t)]);
axis off
hold off;
axis tight;
set(gca,'LooseInset',get(gca,'TightInset'))
end