function obj = compactSpatial(obj)
    for m=1:size(obj.A, 2)
        ai = CNMFE_reshape(obj, obj.A(:, m), 2);
        ai = circular_constraints(ai);
        obj.A(:, m) = ai(:);
    end
end