function positive_ind = temporal_activity_filtering(trace_mat, method, trace_keep_mode)
%% this function is used for assessing the qulity of traces. Two option would be
%  1. use SVM (90% acc and 75% precision in ~4000 traces test dataset)
%  2. use a convolutional neural network.
%  Input
%  trace_mat: matrix which records the temporal traces from captured data
%  last update: 10/26/2020. YZ

% interpolate data
trace_length = 6000;
% check frame step
trace_mat = interp_trace(trace_mat, trace_length);

% make sure the trace mat has specific length
if size(trace_mat, 2) >= trace_length
    trace_mat = trace_mat(:, 1 : trace_length);
else
    trace_mat_new = zeros(size(trace_mat, 1), trace_length);
    trace_mat_new(:, 1 : size(trace_mat, 2)) = trace_mat;
    trace_mat = trace_mat_new;
end

if strcmp(method, 'svm') && strcmp(trace_keep_mode, 'conservative')
    load(sprintf('utility%ssvm_model_sample_6000_conservative.mat', filesep)); 
    XTest_mat_raw = prepare_data_svm(trace_mat, trace_length );
    [YPred_raw, Testscores_raw] = predict(Mdl, XTest_mat_raw);
    positive_ind = find(YPred_raw == categorical(2));
    
elseif strcmp(method, 'cnn') && strcmp(trace_keep_mode, 'conservative')
    load(sprintf('utility%scnn_model_sample_6000_conservative.mat', filesep));    
    % check path length
    miniBatchSize = 50;
    XTest_mat_raw = prepare_data_cnn(trace_mat, trace_length );
    % reshape
    [YPred_raw, Testscores_raw] = classify(net,XTest_mat_raw, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest');
    positive_ind = find(YPred_raw == categorical(2));
elseif strcmp(method, 'svm') && strcmp(trace_keep_mode, 'sensitive')
    load(sprintf('utility%ssvm_model_sample_6000_sensitive.mat', filesep)); 
    XTest_mat_raw = prepare_data_svm(trace_mat, trace_length );
    [YPred_raw, Testscores_raw] = predict(Mdl, XTest_mat_raw);
    positive_ind = find(YPred_raw == categorical(2));
    
elseif strcmp(method, 'cnn') && strcmp(trace_keep_mode, 'sensitive')
    load(sprintf('utility%scnn_model_sample_6000_sensitive.mat', filesep));    
    % check path length
    miniBatchSize = 50;
    XTest_mat_raw = prepare_data_cnn(trace_mat, trace_length );
    % reshape
    [YPred_raw, Testscores_raw] = classify(net,XTest_mat_raw, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest');
    positive_ind = find(YPred_raw == categorical(2));    
else
    error('select svm or cnn as the method!')
end


end

function XTest_mat = prepare_data_cnn(data, datalength)
    assert(datalength == size(data, 2))
    for i = 1 : size(data, 1)
        buf1 = data(i, :);
        % calculate zscore
        buf1 = zscore(buf1, 0, 2);            
        XTest_mat( :, :, 1, i)= buf1(:);   
    end

end
function XTest_mat = prepare_data_svm(data, datalength)
    assert(datalength == size(data, 2))
    for i = 1 : size(data, 1)
        buf1 = data(i, :);
        % calculate zscore
        buf1 = zscore(buf1, 0, 2);            
        XTest_mat(i, :)= buf1(:);   
    end

end

function trace_mat_new = interp_trace(trace_mat, new_length)
    trace_mat_new = imresize(trace_mat, [size(trace_mat, 1), new_length], 'nearest');
end