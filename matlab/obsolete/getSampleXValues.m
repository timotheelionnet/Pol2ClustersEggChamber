function [xSampleVals, xSampleIDs] = getSampleXValues(conditionNames,nSamples,varargin)
    % generates a series of X coordinates for each sample, separated by
    % condition xSampleVals(i,j) is the x coordinate for condition i, sample j 
    % and matching tick names xSampleIDs{i,j}
    
    % conditionNames: 1 x n cell array containing the names of the various
        % conditions.
    % nSamples 1 x n array containing the number of samples per each condition.
    % optional argument spacingUnit is the unit spacing between two samples - its value is irrelevant 
        % unless you want to use the X values in a plot where other things are present and need 
        % to coordinate the X values. default if not populated: 1.
    
    % how much separation there is between conditions as a function of the separation between samples.
    condSeparator = 1.5; 
    
    if numel(varargin) == 0
        spacingUnit = 1;
    else
        spacingUnit = varargin{1};
    end
    
    nConditions = numel(conditionNames);
    
    xSampleVals = zeros(nConditions,max(nSamples));
    xSampleIDs = cell(nConditions,max(nSamples));
    curX = 0;
    for i=1:nConditions
        if i>1
            curX = curX + condSeparator*spacingUnit;
        end
        
        xSampleVals(i,1:nSamples(i)) = curX + (1:nSamples(i))*spacingUnit;
        curX = curX + nSamples(i)*spacingUnit;
    
        for j=1:nSamples(i)
            xSampleIDs{i,j} = [conditionNames{i},'_sample',num2str(j)];
        end
    end


end



