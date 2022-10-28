function t2 = backgroundCorrectNucIntensity(t)
t2 = t;
valueToSubtract = 'Median';
rawSuffix = '_raw';
backgroundIntensityPrefix = {'wholeImg','eggChamber'};

% these variables will get the bg intensity subtracted from them
varsToBeSubtracted = {...
    'Mean',...
    'Max',...
    'Min',...
    'Median',...
    'Mode',...
    };

[nChannels,c] = findChannelNumberFromNucTable(t);

% collect the background columns to subtract
varToSubtract = cell(nChannels,numel(backgroundIntensityPrefix));
for i=1:nChannels
    for j=1:numel(backgroundIntensityPrefix)
        
        varToSubtract{i,j} = [backgroundIntensityPrefix{j},...
            'C',num2str(c(i)),'_',valueToSubtract,rawSuffix];
    end
end

% keep only in t the variables that need subtracting
varList = t.Properties.VariableNames;
idx = [];
for i=1:numel(varsToBeSubtracted)
    curIdx = find(cell2mat( cellfun( @contains, varList,...
        repmat(varsToBeSubtracted(i), size(varList)),...
        'UniformOutput',0) ));
    idx = [idx,curIdx];
end

% remove variables relative to background regions from the 'to be
% subtracted' list
for i=1:numel(backgroundIntensityPrefix)
    curIdx = find(cell2mat( cellfun( @contains, varList,...
        repmat(backgroundIntensityPrefix(i), size(varList)),...
        'UniformOutput',0) ));
    idx = setdiff(idx,curIdx);
end

varList = t.Properties.VariableNames(idx);

% loop through channels and background correct
for i=1:nChannels

    % find variables in t that are relative to the current color channel
    idxC = cell2mat( cellfun( @contains, varList,...
        repmat({['C',num2str(c(i)),'_']}, size(varList)),...
        'UniformOutput',0) );

    curVarList = varList(idxC);
    for j=1:numel(backgroundIntensityPrefix)
        for k =1:numel(curVarList)
            newVarName = strrep(curVarList{k},rawSuffix,['_',backgroundIntensityPrefix{j},'Corr']);
            newVar = t.(curVarList{k}) - t.(varToSubtract{i,j});
            t2 = addvars(t2,newVar,'NewVariableNames',newVarName);
        end
    end
end
