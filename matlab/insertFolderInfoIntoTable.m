function ng = insertFolderInfoIntoTable(ng,idx,folderName,curFolder)
% adds to a nucleus geometry table ng:
    % - a first column with a constant index value (idx) for
        % the current eggChamber (variable 'eggChamberID')
    % - two columns at the end holding the rootFolder where all the 
        % eggchambers subfolders are located (folderName; variable 'rootFolder') and the 
        % current eggChamber sufolder (curFolder; variable 'eggChamberFolder')


% add eggChamber ID in first column
idxVarName = 'eggChamberID';
idxVar = repmat(idx,size(ng,1),1);

ng = addvars(ng,idxVar,'NewVariableNames', idxVarName,'Before',1);

% add root folder in last column
rootVarName = 'rootFolder';
rootVar = repmat(string(folderName),size(ng,1),1);

ng = addvars(ng,rootVar,'NewVariableNames', rootVarName);

% add current folder in last column
curVarName = 'eggChamberFolder';
curVar = repmat(string(curFolder),size(ng,1),1);

ng = addvars(ng,curVar,'NewVariableNames', curVarName);

end