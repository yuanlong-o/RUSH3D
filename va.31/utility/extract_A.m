function A_cell = extract_A(A_mat, movie_size)
A_cell = cell(size(A_mat, 1), movie_size(4));

for i = 1 : size(A_mat, 1) % number of component
    for j = 1 : movie_size(4) % number of view
        curr_A = A_mat((j - 1) * movie_size(1) * movie_size(2) + 1 : j * movie_size(1) * movie_size(2), i);
        A_cell{i, j} = sparse(reshape(curr_A, movie_size(1), movie_size(2)));
    end
end
end