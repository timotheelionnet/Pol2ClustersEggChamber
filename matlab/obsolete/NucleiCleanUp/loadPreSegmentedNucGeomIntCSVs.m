function t = loadPreSegmentedNucGeomIntCSVs(csvFolder)

%% load results into one big table

% collect subfolder names
sd = dir(csvFolder);
ctr = 0;
subDirList = {};
for i=1:numel(sd)
    if isfolder(fullfile(csvFolder,sd(i).name)) ... 
            && ~strcmp(sd(i).name,'.') && ~strcmp(sd(i).name,'..')
        ctr = ctr+1;
        subDirList{ctr,1} = fullfile(csvFolder,sd(i).name);
    end
end

% turn off the warning the Matlab modifies headers that contain
% unacceptable characters (like '.') to decluter the command line.
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% collect files
t = [];
for i=1:numel(subDirList)
    fl = dir(subDirList{i,1});
    for j=1:numel(fl)
        [~,f,e] = fileparts(fl(j).name);
        if strcmp(e,'.csv') && contains(f,'Geom')
            curT = readtable(fullfile(subDirList{i,1},fl(j).name));
            % add conditionID column i.e. index corresponding to the
            % subfolder it came from 
            cond = i*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','conditionID'); 

            % add eggChamberID column i.e. index corresponding to the
            % file it came from 
            cond = j*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','eggChamberID');
            
            %open matching intensity file
            curT2 = readtable(fullfile(subDirList{i,1},strrep(fl(j).name,'Geom','Int')));
            
            curT2 = removevars(curT2,...
                intersect(curT.Properties.VariableNames, curT2.Properties.VariableNames));
            % append to global table
            t = [t;[curT,curT2]];
        end
    end
end


end