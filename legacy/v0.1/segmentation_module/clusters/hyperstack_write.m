function hyperstack_write(file, hyperstack)
% Writes an (up to) 5D stack into a (XYCZT) HyperStack Tiff file
% readable by ImageJ
%
% hyperstack must have class single

% simple checks
assert(nargin == 2, 'Not enough arguments');
% assert(isa(hyperstack, 'single'), 'hyperstack must be single');
if isa(hyperstack, 'single')
    ts.BitsPerSample = 32;
    ts.SampleFormat = Tiff.SampleFormat.IEEEFP;
elseif isa(hyperstack, 'uint8')
	ts.SampleFormat = Tiff.SampleFormat.UInt;
    ts.BitsPerSample = 8;
end

% get all five dimensions
d = zeros(5, 1);
for i = 1 : 5
    d(i) = size(hyperstack, i);
end

% assemble image description
s = sprintf('ImageJ=1.51\nnimages=%d\nchannels=%d\nslices=%d\nframes=%d\nhyperstack=true\nmode=color\nloop=false\nmin=%.1f\nmax=%.1f\n', ...
    prod(d(3:5)), d(3), d(4), d(5), floor(min(hyperstack(:))*10)/10, ceil(max(hyperstack(:))*10)/10);

% open tif file for writing and set file tags
t = Tiff(file, 'w');

ts.ImageLength = d(1);
ts.ImageWidth = d(2);
ts.Photometric = Tiff.Photometric.MinIsBlack;
ts.Compression = Tiff.Compression.None;

ts.SamplesPerPixel = 1;

ts.RowsPerStrip = 5;
ts.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
ts.Software = 'MATLAB';
ts.ImageDescription = s;

% loop over dimensions 3, 4, and 5
for k = 1 : d(5)
    for j = 1 : d(4)
        for i = 1 : d(3)

            frame = hyperstack(:, :, i, j, k);
            t.setTag(ts)            
            t.write(frame);
            t.writeDirectory();
        end
    end
end

% close tif file
t.close();

end