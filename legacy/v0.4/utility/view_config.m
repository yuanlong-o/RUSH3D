function view_array = view_config(preprocess_param)

Nnum = preprocess_param.Nnum;
view_range = preprocess_param.view_range;

[Sx,Sy] = meshgrid(-fix(Nnum/2):fix(Nnum/2),-fix(Nnum/2):fix(Nnum/2));
view_mask = Sx.^2 + Sy.^2 <= view_range^2;
view_array{sum(view_mask(:))} = [];
count = 1;
for u = 1: Nnum
    for v = 1:Nnum
        if view_mask(u,v) == 1
            view_array{count} = [u,v];
            count = count+1;
        end
    end
end
end
