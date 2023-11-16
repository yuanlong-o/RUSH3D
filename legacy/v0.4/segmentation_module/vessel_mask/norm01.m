function normed = norm01(val)
    normed = (val - min(val(:))) / (max(val(:)) - min(val(:)));
end

