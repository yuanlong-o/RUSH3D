function [i_index, j_index] = gen_spiral(Nnum)

if mode(Nnum) == 0
    loop_n = Nnum / 2;
else
    loop_n = (Nnum + 1) / 2;
end
i_index = zeros(1, Nnum * Nnum);
j_index = zeros(1, Nnum * Nnum);

start = 1;
for k = 1 :loop_n
    current_size = Nnum - 2 * (k - 1);
    bias = k - 1;
    if current_size > 1
        i_new = [1 : current_size , current_size * ones(1, current_size - 1), ...
            current_size-1 : -1 : 1, ones(1, current_size - 2)];
        
        j_new = [ones(1, current_size), 2 : current_size, ...
            current_size * ones(1, current_size - 1), current_size - 1 : -1 : 2];
        
        i_index(start : start + 4 * (current_size - 1) - 1) = i_new + bias;
        j_index(start : start + 4 * (current_size - 1) - 1) = j_new + bias;
        start = start + 4 * (current_size - 1);
    else
        i_new = 1;
        j_new = 1;
        i_index(start) = i_new + bias;
        j_index(start) = j_new + bias;        
    end
   
end