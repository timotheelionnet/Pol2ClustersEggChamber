function nr = addImgCorrIntToTable(nr,rawInt,bgInt,c,corrMode,prefix,suffix,insertPosition,addBackgroundValue)
% corrects raw intensity values in rawInt table by subtracting background
% intensity from bgInt, then adds the variables to table nr.
% c is the name of the channel (used to rename variables)
% set corrMode to 'mean' or 'median' to choose which intensity metric to
% subtract.
% prefix is a string to attach at the beginning of the variable name
% suffix is a string to attach at the end of the variable name
% insertPosition is the name of the variable after which to insert the
% first set of corrected variables
% addBackground value is a 0/1 flag to append or not the background value used for
% correction

if size(bgInt,1) == 0
    disp('cannot correct Intensity values for background because background table is empty.');
    return
end

% these variables are inserted in the nr table after the Volume variable
varsToAddFirst = {...
    'Mean',...
    'StdDev',...
    'Max',...
    'Min',...
    'Median',...
    'Mode',...
    'Skewness',...
    'Kurtosis',...
    };

% these variables are inserted in the nr table after the InscrBall_Radius variable
varsToAddLater = {...
    'NumberOfVoxels',...
    'Volume',...
    'NeighborsMean',...
    'NeighborsStdDev',...
    'NeighborsMax',...
    'NeighborsMin',...
    'NeighborsMedian',...
    'NeighborsMode',...
    'NeighborsSkewness',...
    'NeighborsKurtosis',...
    };

% these variables will get the bg intensity subtracted from them
varsToBeSubtracted = {...
    'Mean',...
    'Max',...
    'Min',...
    'Median',...
    'Mode',...
    'NeighborsMean',...
    'NeighborsMax',...
    'NeighborsMin',...
    'NeighborsMedian',...
    'NeighborsMode',...
    };

% remove variables that aren't present in the intensity table
varsToAddFirst = intersect(varsToAddFirst,rawInt.Properties.VariableNames,'stable');
varsToAddLater = intersect(varsToAddLater,rawInt.Properties.VariableNames,'stable');
varsToBeSubtracted = intersect(varsToBeSubtracted,rawInt.Properties.VariableNames,'stable');

% remove variables that don't change upon background subtraction
varsToAddFirst = intersect(varsToAddFirst,varsToBeSubtracted,'stable');
varsToAddLater = intersect(varsToAddLater,varsToBeSubtracted,'stable');

% extract the intensity value to subtract
if strcmp(corrMode,'mean')
    bgVal = bgInt.Mean(1);
else
    bgVal = bgInt.Median(1);
end

% subtract bg intensity to relevant variables in rawInt
for i=1:numel(varsToBeSubtracted)
    rawInt.(varsToBeSubtracted{i}) = rawInt.(varsToBeSubtracted{i}) - bgVal;
end

% change variable names in rawInt
newVarsToAddFirst = varsToAddFirst;
for i=1:numel(varsToAddFirst)
    newVarsToAddFirst{i} = [prefix,'C',num2str(c),'_',newVarsToAddFirst{i},'_',suffix];
end

newVarsToAddLater = varsToAddLater;
for i=1:numel(varsToAddLater)
    newVarsToAddLater{i} = [prefix,'C',num2str(c),'_',newVarsToAddLater{i},'_',suffix];
end

rawInt = renamevars(rawInt,varsToAddFirst,newVarsToAddFirst);
rawInt = renamevars(rawInt,varsToAddLater,newVarsToAddLater);

% insert corrected variables into nr
for i=numel(newVarsToAddFirst):-1:1  % note reverse order because we are adding sequentially to the same position in the table
    nr = addvars(nr,rawInt.(newVarsToAddFirst{i}),...
        'NewVariableNames', newVarsToAddFirst{i},...
        'After',insertPosition);
end

for i=1:numel(newVarsToAddLater)  % append at the end
    nr = addvars(nr,rawInt.(newVarsToAddLater{i}),...
        'NewVariableNames', newVarsToAddLater{i});
end

if addBackgroundValue
    if contains(suffix,'wholeImg')
        bgValPrefix = 'wholeImg';
    elseif contains(suffix,'eggChamber')
        bgValPrefix = 'eggChamber';
    else
        bgValPrefix = ['bgVal',suffix];
    end

    % add variable holding the background intensity value
    bgVarName = [bgValPrefix,'C',num2str(c),'_',corrMode];
    bgVar = repmat(bgVal,size(nr,1),1);
    if ~ismember(bgVarName,nr.Properties.VariableNames)
        nr = addvars(nr,bgVar,'NewVariableNames', bgVarName);
    else
        disp('bg Intensity already present in table, skipping.');
    end
end
end