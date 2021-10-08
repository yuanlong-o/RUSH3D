function [spiral_matrix,seq,seq_ind] = Spiral_circle(side_length,view_range)

[Sx,Sy] = meshgrid(-fix(side_length/2):fix(side_length/2),-fix(side_length/2):fix(side_length/2));
sc_matrix = Sx.^2 + Sy.^2 <= view_range^2;
sc_matrix = double(sc_matrix);
spiral_matrix = zeros(side_length, side_length);
direction = 0; 

init_x = 1; 
init_y = 1;

spiral_matrix(init_x, init_y) = 1;

count = 1;

current_x = init_x;
current_y = init_y;
if sc_matrix(current_x, current_y) ~= 0
    sc_matrix(current_x, current_y) = count;
    count = count + 1;
end

for i = 2 : side_length^2
    
    if direction == 0
        
        current_y = current_y + 1;
        spiral_matrix(current_x, current_y) = i;
        if sc_matrix(current_x, current_y) ~= 0
            sc_matrix(current_x, current_y) = count;
            count = count + 1;
        end
        if current_y + 1 > side_length || spiral_matrix(current_x, current_y + 1) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    elseif direction == 1
        
        current_x = current_x + 1;
        spiral_matrix(current_x, current_y) = i;
        if sc_matrix(current_x, current_y) ~= 0
            sc_matrix(current_x, current_y) = count;
            count = count + 1;
        end
        if current_x + 1 > side_length || spiral_matrix(current_x + 1, current_y) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    elseif direction == 2
        
        current_y = current_y - 1;
        spiral_matrix(current_x, current_y) = i;
        if sc_matrix(current_x, current_y) ~= 0
            sc_matrix(current_x, current_y) = count;
            count = count + 1;
        end
        if current_y - 1 < 1 || spiral_matrix(current_x, current_y - 1) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    else 
        
        current_x = current_x - 1;
        spiral_matrix(current_x, current_y) = i;
        if sc_matrix(current_x, current_y) ~= 0
            sc_matrix(current_x, current_y) = count;
            count = count + 1;
        end
        if (current_x - 1 < 1) || (spiral_matrix(current_x - 1, current_y) ~= 0)
            direction = mod(direction + 1, 4);
        end 
    end
end

seq = zeros(1,count-1);
seq_ind = zeros(1,count-1);
count = 1;
for u = 1: side_length
    for v = 1:side_length
        if sc_matrix(u,v) ~= 0
            seq(count) = sc_matrix(u,v);
            count = count+1;
        end
    end
end

for seq_id = 1: length(seq)
    seq_ind(seq_id) = find(seq_id == seq);
end
end
