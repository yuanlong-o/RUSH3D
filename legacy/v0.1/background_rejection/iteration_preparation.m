function processed_video = iteration_preparation(realigndata_name_perfix, valid_frame, specify_wigner)
%  load wigners, and subtract marginal view information from central view
%  output is the processed video.
%  last update: 5/23/2021. YZ


num_wigner = size(specify_wigner, 1);

buf_wigner_video = loadtiff(sprintf('%s_No%d.tif', realigndata_name_perfix, 1));
[size_h, size_w, ~] = size(buf_wigner_video);

processed_video = zeros(size_h, size_w, valid_frame, num_wigner, 'single');
for i = 1 : valid_frame
    if mod(i, 10) == 0
       fprintf('%d in %d frames collected \n', i, valid_frame) 
    end
    curr_video = loadtiff(sprintf('%s_No%d.tif', realigndata_name_perfix, i - 1));
    curr_video = single(curr_video) / 65535;
    for j = 1 : num_wigner
        processed_video(:, :, i, j) = curr_video(:, :, specify_wigner(j));
    end
end


end