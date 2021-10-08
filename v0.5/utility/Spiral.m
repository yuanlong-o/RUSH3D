function spiral_matrix = Spiral(side_length)
spiral_matrix = zeros(side_length, side_length);
direction = 0; 

init_x = 1; 
init_y = 1;

spiral_matrix(init_x, init_y) = 1;

current_x = init_x;
current_y = init_y;

for i = 2 : side_length^2
    
    if direction == 0
        
        current_y = current_y + 1;
        spiral_matrix(current_x, current_y) = i;
        if current_y + 1 > side_length || spiral_matrix(current_x, current_y + 1) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    elseif direction == 1
        
        current_x = current_x + 1;
        spiral_matrix(current_x, current_y) = i;
        if current_x + 1 > side_length || spiral_matrix(current_x + 1, current_y) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    elseif direction == 2
        
        current_y = current_y - 1;
        spiral_matrix(current_x, current_y) = i;
        if current_y - 1 < 1 || spiral_matrix(current_x, current_y - 1) ~= 0
            direction = mod(direction + 1, 4);
        end
        
    else 
        
        current_x = current_x - 1;
        spiral_matrix(current_x, current_y) = i;
        if (current_x - 1 < 1) || (spiral_matrix(current_x - 1, current_y) ~= 0)
            direction = mod(direction + 1, 4);
        end
        
    end
end
end
