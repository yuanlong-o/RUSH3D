%% utility function
function out_mat_A = cell_process_A(A, movie_size) % would be slow since reach the limit
out_mat_A = zeros(movie_size(1) * movie_size(2) * movie_size(4), size(A, 1));
%out_mat_A = zeros(movie_size(1) * movie_size(2) * movie_size(4), size(A, 2));
for i = 1 : size(A, 1) % number of component
    %         i
    curr_A =  cell2mat(A(i, :));
    out_mat_A(:, i) = curr_A(:);
end
out_mat_A  = sparse(out_mat_A );
end