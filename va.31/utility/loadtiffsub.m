function I_sub_stack = loadtiffsub(filename, tl, br)
Iinfo = imfinfo(filename);
h = br(1) - tl(1) + 1;
w = br(2) - tl(2) + 1;
I_sub_stack = zeros(h,w,length(Iinfo));
for t = 1: length(Iinfo)
    I_sub_stack(:,:,t) = double(imread(filename, t, 'PixelRegion',{[tl(1),br(1)],[tl(2),br(2)]}));
end
end