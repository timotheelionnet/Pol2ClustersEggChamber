function clusterTable = loadClusterGeom(folderName,curFolder,eggChamberID,nucleusID,uniqueID)
    
    csvFolder = fullfile(fullfile(folderName,curFolder),'csvFIJI');
    
    % load first channel to extract volume 
    fileName = ['nuc',num2str(nucleusID),'_clustGeom.csv'];
    clusterTable = readtable(fullfile(csvFolder,fileName));

    if size(clusterTable,1) == 0
        clusterTable = [];
        return
    end
    
    % add modifier to coordinates to differentiate local vs global
    % coordinates
    oldVars = {'Box.X.Min','Box.X.Max','Box.Y.Min',...
        'Box.Y.Max','Box.Z.Min','Box.Z.Max',...
        'Centroid.X','Centroid.Y','Centroid.Z',...
        'Elli.Center.X','Elli.Center.Y','Elli.Center.Z',...
        'InscrBall.Center.X','InscrBall.Center.Y','InscrBall.Center.Z'};
    oldVars = cellfun(@strrep,oldVars,...
        repmat({'.'},size(oldVars)),repmat({'_'},size(oldVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,oldVars,...
        repmat({'X'},size(oldVars)),repmat({'localX'},size(oldVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,newVars,...
        repmat({'Y'},size(newVars)),repmat({'localY'},size(newVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,newVars,...
        repmat({'Z'},size(newVars)),repmat({'localZ'},size(newVars)), ...
        'UniformOutput',0);

    for i=1:numel(clusterTable.Properties.VariableNames)
        if ismember(clusterTable.Properties.VariableNames{i},oldVars)
            idx = ismember(oldVars,clusterTable.Properties.VariableNames{i});
            clusterTable.Properties.VariableNames{i} = ...
                ['clust',newVars{idx}];
        else
            clusterTable.Properties.VariableNames{i} = ...
                ['clust',clusterTable.Properties.VariableNames{i}];
        end
    end

    % add  nucleus ID in first column (will become third column after
    % we add the other IDs)
    idxVarName = 'nucID';
    idxVar = repmat(nucleusID,size(clusterTable,1),1);
    
    clusterTable = addvars(clusterTable,idxVar,'NewVariableNames', idxVarName,'Before',1);

    % add  eggChamber ID in first column
    idxVarName = 'eggChamberID';
    idxVar = repmat(eggChamberID,size(clusterTable,1),1);
    
    clusterTable = addvars(clusterTable,idxVar,'NewVariableNames', idxVarName,'Before',1);
    
    % add  unqiue ID in first column
    idxVarName = 'nucUniqueID';
    idxVar = repmat(uniqueID,size(clusterTable,1),1);
    
    clusterTable = addvars(clusterTable,idxVar,'NewVariableNames', idxVarName,'Before',1);
    

    % add root folder in last column
    rootVarName = 'rootFolder';
    rootVar = repmat(string(folderName),size(clusterTable,1),1);
    
    clusterTable = addvars(clusterTable,rootVar,'NewVariableNames', rootVarName);
    
    % add current folder in last column
    curVarName = 'eggChamberFolder';
    curVar = repmat(string(curFolder),size(clusterTable,1),1);
    
    clusterTable = addvars(clusterTable,curVar,'NewVariableNames', curVarName);


end