function [tuple, hashkeylist] = computeTuples(inputImage,operator)

if ~iscell(operator)
    operatorlist{1} = operator;
else
    operatorlist = operator;
end

set = cell(0);
index = 1;

for op = 1:length(operatorlist)
    params = operatorlist{op}.params;
    reflectOperator.tuples = operatorlist{op}.tuples;
                 
    for reflection = 1:2^params.invariance.reflection    
        if reflection == 2
            if strcmp(params.inputfilter.name,'Gabor')
                if params.inputfilter.symmetric == 1                
                    reflectOperator.tuples(2,:) = mod(pi - reflectOperator.tuples(2,:),pi);                
                else
                    reflectOperator.tuples(2,:) = mod(pi - reflectOperator.tuples(2,:),2*pi);
                end
            end
            reflectOperator.tuples(4,:) = mod(pi - reflectOperator.tuples(4,:),2*pi);
        end

        rotateOperator = reflectOperator;
        for psiindex = 1:length(params.invariance.rotation.psilist)
            if strcmp(params.inputfilter.name,'Gabor')
                if params.inputfilter.symmetric == 1
                    rotateOperator.tuples(2,:) = mod(reflectOperator.tuples(2,:) + params.invariance.rotation.psilist(psiindex),pi);
                else
                    rotateOperator.tuples(2,:) = mod(reflectOperator.tuples(2,:) + params.invariance.rotation.psilist(psiindex),2*pi);
                end
            end
            rotateOperator.tuples(4,:) = mod(reflectOperator.tuples(4,:) + params.invariance.rotation.psilist(psiindex),2*pi);            

            scaleOperator = rotateOperator;
            for upsilonindex = 1:length(params.invariance.scale.upsilonlist)
                % Scale the values of parameters lambda and rho of every tuple by a given upsilon value
                
                if strcmp(params.inputfilter.name,'Gabor')
                    scaleOperator.tuples(1,:) = rotateOperator.tuples(1,:) * params.invariance.scale.upsilonlist(upsilonindex);   
                elseif strcmp(params.inputfilter.name,'DoG')
                    scaleOperator.tuples(2,:) = rotateOperator.tuples(2,:) * params.invariance.scale.upsilonlist(upsilonindex);   
                end
                
                scaleOperator.tuples(3,:) = rotateOperator.tuples(3,:) * params.invariance.scale.upsilonlist(upsilonindex);

                if strcmp(params.COSFIRE.outputfunction,'weightedgeometricmean')
                    tupleweightsigma = sqrt(-max(scaleOperator.tuples(3,:))^2/(2*log(params.COSFIRE.mintupleweight)));
                    scaleOperator.tuples(5,:) = exp(-scaleOperator.tuples(3,:).^2./(2*tupleweightsigma*tupleweightsigma));
                else
                    scaleOperator.tuples(5,:) = ones(size(scaleOperator.tuples(3,:))); 
                end
                scaleOperator.tuples(6,:) = params.COSFIRE.sigma0; 
                scaleOperator.tuples(7,:) = params.COSFIRE.alpha; 
                set{index} = scaleOperator.tuples;
                index = index + 1;
            end    
        end
    end
end

S = cell2mat(set);
S = round(S * 10000) / 10000;
inputfilterparams = unique(S(1:2,:)','rows');
tupleparams = unique(S([1:3,5,6,7],:)','rows');
sz = size(inputImage);

switch (params.ht)
    case 0        
        inputfilterresponse = zeros(sz(1),sz(2),size(inputfilterparams,1));
        for i = 1:size(inputfilterparams,1)
            inputfilterresponse(:,:,i) = getGaborResponse(inputImage,params,inputfilterparams(i,1),inputfilterparams(i,2));    
        end
        inputfilterresponse(inputfilterresponse < params.COSFIRE.t1*max(inputfilterresponse(:))) = 0;
        
        response = cell(1,size(tupleparams,1));
        for i = 1:size(tupleparams,1)
            index = ismember(inputfilterparams,tupleparams(i,1:2),'rows');
            ifresp = inputfilterresponse(:,:,index);
            
            % Compute the sigma of the 2D Gaussian function that will be
            % used to blur the corresponding output.
            sigma = (params.COSFIRE.sigma0 + (params.COSFIRE.alpha*tupleparams(i,3)));                

            if any(ifresp(:))
                if strcmp(params.COSFIRE.blurringfunction,'max')
                    response{i} = blurshift(ifresp,sigma,0,0);
                elseif strcmp(params.COSFIRE.blurringfunction,'sum')
                    blurfunction = fspecial('gaussian',round([sigma sigma].*6),sigma);
                    response{i} = conv2(ifresp,blurfunction,'same');
                    response{i} = circshift(response{i},[0,0]);
                end
            else
                response{i} = zeros(sz);
            end            

            if strcmp(params.COSFIRE.outputfunction,'weightedgeometricmean')                     
                weight = exp(-tupleparams(i,3)^2/(2*params.COSFIRE.weightingsigma*params.COSFIRE.weightingsigma)); 
                response{i} = response{i} .^ weight;                    
            end
        end
        tuple.response = response;
        tuple.params = tupleparams;          
    case 1
        inputfilterresponse = zeros(size(inputImage,1),size(inputImage,2),size(inputfilterparams,1));
        ninputfilterparams = size(inputfilterparams,1);
        hashkeylist = cell(1,ninputfilterparams);
        hashvaluelist = cell(1,ninputfilterparams);
        
        % Compute the responses of the bank of input filters
        for i = 1:ninputfilterparams
            hashkeylist{i} = getHashkey(inputfilterparams(i,:));
            hashvaluelist{i} = i;
            
            if strcmp(params.inputfilter.name,'Gabor')
                inputfilterresponse(:,:,i) = getGaborResponse(inputImage,params.inputfilter.Gabor,inputfilterparams(i,1),inputfilterparams(i,2));    
            elseif strcmp(params.inputfilter.name,'DoG')
                inputfilterresponse(:,:,i) = getDoG(inputImage,inputfilterparams(i,2), inputfilterparams(i,1), params.inputfilter.DoG.sigmaratio, 0, params.inputfilter.DoG.halfwaverect);
            end
        end
        DOG = inputfilterresponse(:,:,i);
        [num idx] = max(DOG(:));
        [y x] = ind2sub(size(DOG), idx);
        
        inputfilterhashtable = containers.Map(hashkeylist,hashvaluelist);
        
        % Threshold the responses by a fraction t1 of the maximum response
        % inputfilterresponse(inputfilterresponse < params.COSFIRE.t1*max(inputfilterresponse(:))) = 0;
        
        ntupleparams = size(tupleparams,1);
        hashkeylist = cell(1,ntupleparams);
        hashvaluelist = cell(1,ntupleparams);

        for i = 1:ntupleparams
            rho    = tupleparams(i,3);
            index = inputfilterhashtable(getHashkey(tupleparams(i,1:2)));
            hashkeylist{i} = getHashkey(tupleparams(i,[1:3,5,6]));    
            ifresp = inputfilterresponse(:,:,index);
            
            if any(ifresp(:))               
                % Compute the sigma of the 2D Gaussian function that will be
                % used to blur the corresponding output.
                %sigma = (params.COSFIRE.sigma0 + (params.COSFIRE.alpha*rho));                
                sigma = (tupleparams(i,5) + (tupleparams(i,6) * rho));                
                %tic;
                if strcmp(params.COSFIRE.blurringfunction,'max')                    
                    hashvaluelist{i} = blurshift(ifresp,sigma,0,0); 
                elseif strcmp(params.COSFIRE.blurringfunction,'sum')
                    blurfunction = fspecial('gaussian',[1 round(sigma.*6)],sigma);
                    hashvaluelist{i} = conv2(blurfunction,blurfunction,ifresp,'same');
                end
                %toc;
                
                BLUR = hashvaluelist{i};
                [num idx] = max(BLUR(:));
                [y x] = ind2sub(size(BLUR), idx);
            else
                hashvaluelist{i} = zeros(sz);
            end

            if strcmp(params.COSFIRE.outputfunction,'weightedgeometricmean')                                     
                hashvaluelist{i} = hashvaluelist{i} .^ tupleparams(i,4);                    
            end    
        end
        [~,srtidx] = sort(cell2mat(hashkeylist));
        tuple.hashtable  = containers.Map(hashkeylist(srtidx),hashvaluelist(srtidx));  
end