function [resp_stack] = BCOSFIRE_lfm(image_stack, filter1, filter2)
% Delineation of blood vessels in retinal images based on combination of BCOSFIRE filters responses.
%
% VERSION 02/09/2016
% CREATED BY: George Azzopardi (1), Nicola Strisciuglio (1,2), Mario Vento (2) and Nicolai Petkov (1)
%             1) University of Groningen, Johann Bernoulli Institute for Mathematics and Computer Science, Intelligent Systems
%             1) University of Salerno, Dept. of Information Eng., Electrical Eng. and Applied Math., MIVIA Lab
%
%   If you use this script please cite the following paper:
%   [1] "George Azzopardi, Nicola Strisciuglio, Mario Vento, Nicolai Petkov,
%   Trainable COSFIRE filters for vessel delineation with application to retinal images,
%   Medical Image Analysis, Volume 19 , Issue 1 , 46 - 57, ISSN 1361-8415,
%   http://dx.doi.org/10.1016/j.media.2014.08.002"
%
%   BCOSFIRE_media15 achieves orientation selectivity by combining the output - at certain
%   positions with respect to the center of the COSFIRE filter - of center-on
%   difference of Gaussians (DoG) functions by a  geometric mean.
%
%   BCOSFIRE_media15 takes as input:
%      image -> RGB retinal fundus image
%      filter1 -> a struct defining the configuration parameters of the
%                 symmetric filter:
%                   Name                Value
%                   'sigma'             Value of the standard deviation of
%                                       the outer Gaussian funciton in the DoG
%                   'len'               The length of the filter support
%                   'sigma0'            Value of the standard deviation of
%                                       the Gaussian weigthing function
%                   'alpha'             Alpha coefficient of the weighting
%                                       function
%      filter2 -> a struct defining the configuration parameters of the
%                 asymmetric filter:
%                   Name                Value
%                   'sigma'             Value of the standard deviation of
%                                       the outer Gaussian funciton in the DoG
%                   'len'               The length of the filter support
%                   'sigma0'            Value of the standard deviation of
%                                       the Gaussian weigthing function
%                   'alpha'             Alpha coefficient of the weighting
%                                       function
%      preprocessthresh -> a threshold value used in the preprocessing step
%                          (must be a float value in [0, 1])
%      thresh -> threshold value used to threshold the final filter
%                response (must be an integer value in [0, 255])
%
%   BCOSFIRE returns:
%      resp         -> response of the combination of a symmetric and an
%                   asymemtric COSFIRE filters
%      oriensmap    -> map of the orientation that gives the strongest
%                   response for each pixel.
%
%   The ouput parameter 'oriensmap' is optional. In case it is not
%   required, the algorithm provides a response image as described in [1].
%   On the contrary, the response image is computed by first summing up the
%   symmetric and asymmetric B-COSFIRE filter responses at each
%   orientation, and then superimposing such responses to achieve rotation
%   invariance. In this case, the orientation map, with information
%   about the orientation that gives the strongest response at every pixel,
%   is provided as output.
%
%   Examples: [resp] = BCOSFIRE(imread('01_test.tif'), f1, f2, 0.5, 38);
%             [resp oriensmap] = BCOSFIRE(imread('01_test.tif'), f1, f2, 0.5, 38);
%
%   The image 01_test.tif is taken from the DRIVE data set, which can be
%   downloaded from: http://www.isi.uu.nl/Research/Databases/DRIVE/

% path(path,'./Gabor/');
% path(path,'./COSFIRE/');
% path(path,'./Preprocessing/');
% path(path,'./Performance/');

%% Model configuration
% Prototype pattern
x = 101; y = 101; % center
line1(:, :) = zeros(201);
line1(:, x) = 1; %prototype line

% Symmetric filter params
symmfilter = cell(1);
symm_params = SystemConfig;
% COSFIRE params
symm_params.inputfilter.DoG.sigmalist = filter1.sigma;
symm_params.COSFIRE.rholist = 0:2:filter1.len;
symm_params.COSFIRE.sigma0 = filter1.sigma0 / 6;
symm_params.COSFIRE.alpha = filter1.alpha / 6;
% Orientations
numoriens = 12;
symm_params.invariance.rotation.psilist = 0:pi/numoriens:pi-pi/numoriens;
% Configuration
symmfilter{1} = configureCOSFIRE(line1, round([y x]), symm_params);
% Show the structure of the COSFIRE filter
% showCOSFIREstructure(symmfilter);

if filter2
    % Asymmetric filter params
    asymmfilter = cell(1);
    asymm_params = SystemConfig;
    % COSFIRE params
    asymm_params.inputfilter.DoG.sigmalist = filter2.sigma;
    asymm_params.COSFIRE.rholist = 0:2:filter2.len;
    asymm_params.COSFIRE.sigma0 = filter2.sigma0 / 6;
    asymm_params.COSFIRE.alpha = filter2.alpha / 6;
    % Orientations
    numoriens = 24;
    asymm_params.invariance.rotation.psilist = 0:2*pi/numoriens:(2*pi)-(2*pi/numoriens);
    % Configuration
    asymmfilter{1} = configureCOSFIRE(line1, round([y x]), asymm_params);
    asymmfilter{1}.tuples(:, asymmfilter{1}.tuples(4,:) > pi) = []; % Deletion of on side of the filter
end
% Show the structure of the COSFIRE operator
% showCOSFIREstructure(asymmfilter);

%% Filtering
resp_stack = zeros(size(image_stack));
% textprogressbar('Applying BCOSFIRE filter ');
for zix = 1:size(image_stack,3)
    image = preprocess_lfm(squeeze(image_stack(:,:,zix)));
    %image = 1 - image;
    
    % Apply the symmetric B-COSFIRE to the input image
    rot1 = applyCOSFIRE(image, symmfilter);
    if filter2
        rot2 = applyCOSFIRE(image, asymmfilter);
    end
    %
    % if nargout == 1
    % The code as presented in the paper
    rot1 = max(rot1{1},[],3);
    if filter2
        rot2 = max(rot2{1},[],3);
        resp = rot1 + rot2;
    else
        resp = rot1;
    end
    % elseif nargout == 2
    %     % Modified code to also give the orientation map as output
    %     for i = 1:size(rot1{1},3)
    %         resp(:,:,i) = rot1{1}(:, :, i) + max(rot2{1}(:,:,i),rot2{1}(:,:,i+12));
    %     end
    %     [resp,oriensmap] = max(resp, [], 3);
    %     oriensmap = symm_params.invariance.rotation.psilist(oriensmap);
    %     oriensmap = oriensmap .* mask;
    % end
    
    resp_stack(:,:,zix) = norm01(resp);
%     textprogressbar(zix/size(image_stack,3) * 100);
end
% textprogressbar(' done');
end

function [img] = preprocess_lfm(img)
mask = ones(size(img));
%img(img < 0.1 * max(img(:))) = 0;
[ignore img] = getBigimg(norm01(img), mask);
img = adapthisteq(img);
end

function [bigimg smallimg] = getBigimg(img,mask)

[sizey, sizex] = size(img);

bigimg = zeros(sizey + 100, sizex + 100);
bigimg(51:(50+sizey), 51:(50+sizex)) = img;

bigmask = false(sizey + 100, sizex + 100);
bigmask(51:(50+sizey), (51:50+sizex)) = mask;

% Creates artificial extension of image.
bigimg = fakepad(bigimg, bigmask, 5, 10);
smallimg = bigimg(51:(50+sizey), 51:(50+sizex));
end
