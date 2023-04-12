function cl = addClusterGlobalCoordinatesToTable(cl,ng,nucIdx)

% isolate the geometry variables that are affected by a change of coordinates
% references
varsToOffset = {'Box.X.Min','Box.X.Max',...
    'Box.Y.Min','Box.Y.Max','Box.Z.Min','Box.Z.Max',...
    'Centroid.X','Centroid.Y','Centroid.Z',...
    'Elli.Center.X','Elli.Center.Y','Elli.Center.Z',...
    'InscrBall.Center.X','InscrBall.Center.Y','InscrBall.Center.Z'};

oldVars = cellfun(@strrep,varsToOffset,...
        repmat({'.'},size(varsToOffset)),repmat({'_'},size(varsToOffset)), ...
        'UniformOutput',0);

oldVars = cellfun(@strrep,oldVars,...
    repmat({'X'},size(oldVars)),repmat({'localX'},size(oldVars)), ...
    'UniformOutput',0);

oldVars = cellfun(@strrep,oldVars,...
    repmat({'Y'},size(oldVars)),repmat({'localY'},size(oldVars)), ...
    'UniformOutput',0);

oldVars = cellfun(@strrep,oldVars,...
    repmat({'Z'},size(oldVars)),repmat({'localZ'},size(oldVars)), ...
    'UniformOutput',0);

for i=1:numel(oldVars)
    oldVars{i} = ['clust',oldVars{i}];
end

% set up a list of new variable names
newVars = cellfun(@strrep,oldVars,...
    repmat({'local'},size(oldVars)),repmat({''},size(oldVars)), ...
    'UniformOutput',0);

% compute the transformation offset from local coordinates (relative to cropped nucleus stack)
% to global coordinates (relative to whole egg chamber)
dX = ng.nucCentroid_X(nucIdx) - ng.nucCentroid_localX(nucIdx);
dY = ng.nucCentroid_Y(nucIdx) - ng.nucCentroid_localY(nucIdx);
dZ = ng.nucCentroid_Z(nucIdx) - ng.nucCentroid_localZ(nucIdx);

for i=1:numel(oldVars)
    idx = ismember(cl.Properties.VariableNames,oldVars{i});
    if(sum(idx)==0)
        disp(['Could not find variable ',oldVars{i},' in current cluster table.']);
    else
        if contains(oldVars{i},'localX')
            offsetVar = cl.(oldVars{i}) + dX; 
        elseif contains(oldVars{i},'localY')
            offsetVar = cl.(oldVars{i}) + dY; 
        elseif contains(oldVars{i},'localZ')
            offsetVar = cl.(oldVars{i}) + dZ; 
        else
            disp(['Could not find which coordinate variable ',oldVars{i},...
                ' corresponds to (X,Y,Z?) in current cluster table.']);
            offsetVar = [];
        end
        cl = addvars(cl,offsetVar,'NewVariableNames',{newVars{i}});
    end
end


end



