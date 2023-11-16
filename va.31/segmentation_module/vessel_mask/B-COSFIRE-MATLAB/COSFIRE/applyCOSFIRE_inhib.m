function [output, rotations] = applyCOSFIRE_inhib(inputImage, operatorlist, inhibfactor, tuple)
% Delineation of elongated patterns in images.
% Application of a B-COSFIRE filter with or without inhibition
%
%   version: 03/07/2017
%   author: Nicola Strisciuglio
%
%
%   If you use this script please cite the following paper:
%   "George Azzopardi, Nicola Strisciuglio, Mario Vento, Nicolai Petkov, 
%   Trainable COSFIRE filters for vessel delineation with application to retinal images, 
%   Medical Image Analysis, Volume 19 , Issue 1 , 46 - 57, ISSN 1361-8415, 
%   http://dx.doi.org/10.1016/j.media.2014.08.002"


    
    if nargin < 4
        tuple = computeTuples(inputImage, operatorlist);
    end
    output = cell(1); % = cell(1,length(operatorlist));

    excitatoryOperator = operatorlist{1};
    if inhibfactor > 0
        inhibitoryOperator = operatorlist{2};
    else
        inhibitoryOperator = [];
    end

    % COSFIRE response computation
    cosfirefun = @(inputImage,operator,tuple) computeCOSFIRE(inputImage, operator, tuple);                                                  
%     scalefun = @(inputImage,operator,tuple) scaleInvariantCOSFIRE(inputImage,excitatoryOperator,inhibitoryOperator,tuple,cosfirefun);
%     rotationfun = @(inputImage,operator,tuple) rotationInvariantCOSFIRE(inputImage,excitatoryOperator,inhibitoryOperator,tuple,scalefun);
%     output{opindex} = reflectionInvariantCOSFIRE(inputImage,excitatoryOperator,inhibitoryOperator,tuple,rotationfun);                     

    % We compute COSFIRE responses only in rotation-tolerant mode
    if nargout == 2
        [output{1}, rotations] = rotationInvariantCOSFIRE(inputImage, excitatoryOperator, inhibitoryOperator, tuple, cosfirefun, inhibfactor);                     
    else
        output{1} = rotationInvariantCOSFIRE(inputImage, excitatoryOperator, inhibitoryOperator, tuple, cosfirefun, inhibfactor);                     
    end
    
    % Suppress values that are less than a fraction t3 of the maximum
    %output{1}(output{1} < operator.params.COSFIRE.t3 * max(output{1}(:))) = 0;
%end

% function [output] = reflectionInvariantCOSFIRE(inputImage,operator,operator2,tuple,funCOSFIRE)
% 
% Apply the given COSFIRE filter 
% output = feval(funCOSFIRE,inputImage,operator,tuple);
% 
% if operator.params.invariance.reflection == 1
%     Apply a COSFIRE filter which is selective for a reflected version about the y-axis 
%     of the pattern of interest 
%     reflectionDetector = operator;
% 
%     if strcmp(operator.params.inputfilter.name,'Gabor')
%         if operator.params.inputfilter.symmetric == 1
%             reflectionDetector.tuples(2,:) = mod(pi - reflectionDetector.tuples(2,:),pi);
%         else
%             reflectionDetector.tuples(2,:) = mod(pi - reflectionDetector.tuples(2,:),2*pi);
%         end
%     end
%     reflectionDetector.tuples(4,:) = mod(pi - reflectionDetector.tuples(4,:),2*pi);
% 
%     reflectionoutput = feval(funCOSFIRE,inputImage,reflectionDetector,tuple);
% 
%     Take the maximum value of the output of the two COSFIRE filters
%     output = max(output, reflectionoutput);
% end

function [output rotations] = rotationInvariantCOSFIRE(inputImage, excitatoryOperator, inhibitoryOperator, tuple, funCOSFIRE, inhibfactor)

    output = zeros(size(inputImage));
    noriens = length(excitatoryOperator.params.invariance.rotation.psilist);
    if nargout == 2
        rotations = zeros([size(inputImage) noriens]);
    end

    rotateExcitatoryDetector = excitatoryOperator;
    rotateInhibitoryDetector = inhibitoryOperator;

    for psiindex = 1:noriens 
        % Shift the values of parameters (theta,rho) of every tuple by a given psi value
        rotateExcitatoryDetector.tuples(4,:) = excitatoryOperator.tuples(4,:) + excitatoryOperator.params.invariance.rotation.psilist(psiindex);            
        % Compute the output of COSFIRE for the given psi value
        rotexcoutput = feval(funCOSFIRE,inputImage,rotateExcitatoryDetector,tuple);    
        
        % Compute inhibitory response (if enabled)
        if inhibfactor > 0
            rotateInhibitoryDetector.tuples(4,:) = inhibitoryOperator.tuples(4,:) + inhibitoryOperator.params.invariance.rotation.psilist(psiindex);            
            rotinhoutput = feval(funCOSFIRE, inputImage, rotateInhibitoryDetector, tuple);    
            rotoutput = rotexcoutput - (inhibfactor .* rotinhoutput);
            rotoutput = (rotoutput > 0) .* rotoutput;
        else
            rotoutput = rotexcoutput;
        end
        
        if nargout == 2
            rotations(:,:, psiindex) = rotoutput;
        end
        
        % Take the maximum over the COSFIRE outputs for all given values of psi
        output = max(rotoutput, output);
    end

% function [output] = scaleInvariantCOSFIRE(inputImage,operator,operator2,tuple,funCOSFIRE)
% 
% output = zeros(size(inputImage));
% scaleDetector = operator;
% 
% for upsilonindex = 1:length(operator.params.invariance.scale.upsilonlist)
%     Scale the values of parameters lambda and rho of every tuple by a given upsilon value
%     
%     if strcmp(operator.params.inputfilter.name,'Gabor')
%         scaleDetector.tuples(1,:) = operator.tuples(1,:) * operator.params.invariance.scale.upsilonlist(upsilonindex);   
%     elseif strcmp(operator.params.inputfilter.name,'DoG')
%         scaleDetector.tuples(2,:) = operator.tuples(2,:) * operator.params.invariance.scale.upsilonlist(upsilonindex);   
%     end
%     scaleDetector.tuples(3,:) = operator.tuples(3,:) * operator.params.invariance.scale.upsilonlist(upsilonindex);
%         
%     Compute the output of COSFIRE for the given psi value
%     scaleoutput = feval(funCOSFIRE,inputImage,scaleDetector,tuple);
%     
%     Take the maximum over the COSFIRE outputs for all given values of upsilon
%     output = max(output,scaleoutput);
% end

function [output] = computeCOSFIRE(inputImage,operator,tuple)       
    operator.tuples = round(operator.tuples * 10000) / 10000;       
    sigma0 = round(operator.params.COSFIRE.sigma0 * 10000) / 10000;       
    alpha = round(operator.params.COSFIRE.alpha * 10000) / 10000;       

    sz = size(inputImage);
    output = ones(sz);
    ntuples = size(operator.tuples,2);
    %tupleoutputs = 0;
    % Loop through all tuples (subunits) of the operator
    outputs = zeros(sz(1), sz(2), ntuples);
    for sindex = 1:ntuples
        % Convert the polar-coordinate shift vector (rho,phi) to image coordinates
        [col row] = pol2cart(operator.tuples(4,sindex),operator.tuples(3,sindex));  
       
        switch (operator.params.ht)
            case 0
                index = ismember(tuple.params,operator.tuples(1:3,sindex)','rows');
                % Bugfix #1: approximation of the position was not correct.
                % fix() was substituted with round()
                tupleoutput = circshift(tuple.response{index},[fix(row), -fix(col)]); 
            case 1
                hashkey = getHashkey([operator.tuples(1:3,sindex)',sigma0,alpha]);
                % Bugfix #1: approximation of the position was not correct.
                % fix() was substituted with round()
                tupleoutput = circshift(tuple.hashtable(hashkey),[fix(row), -fix(col)]); 
        end

        outputs(:, :, sindex) = tupleoutput; % intermediate responses
        output = output .* tupleoutput; % construction of final response
        
        if ~any(output(:))
            output = zeros(sz);
            return;
        end    
    end

    if strcmp(operator.params.COSFIRE.outputfunction, 'weightedgeometricmean')
        % Compute the COSFIRE output using weighted geometric mean
        tupleweightsigma = sqrt(-max(operator.tuples(3,:))^2/(2*log(operator.params.COSFIRE.mintupleweight)));
        tupleweight = exp(-(operator.tuples(3,:).^2)./(2*tupleweightsigma*tupleweightsigma));    
        output = output .^ (1/sum(tupleweight));    
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'min')
        % Compute the COSFIRE output using min operation
        output = min(outputs, [], 3);
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'geometricmean')
        % Compute the COSFIRE output using geometric mean
        m = output .^ (1/ntuples);
        output = m;
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'arithmeticmean')
        output = mean(outputs, 3);
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'harmonicmean')
        output = harmmean(outputs, 3);
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'trimmedmean')
        mn = min(outputs,[],3);
        ind = mn == 0;
        output = trimmean(outputs, 80, 'weighted', 3);
        output(ind) = 0;
    elseif strcmp(operator.params.COSFIRE.outputfunction, 'median')
        mn = min(outputs,[],3);
        ind = mn == 0;
        output = median(outputs, 3);
        output(ind) = 0;
    else
        % Other combination functions can be added here.
    end
    outputs = [];