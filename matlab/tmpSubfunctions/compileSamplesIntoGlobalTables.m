folderName = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/20221003/res';
% folderName = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/20221003/res';
% collects all nuclei-relative results in one table 'allNuclei.txt"
% collects all cluster-relative results in one table 'allCLusters.txt"


% collect list of subfolders
d = dir(folderName);
dFolders = d([d(:).isdir]);
dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));

% setting that tells the script to use the median of the image / or the egg
% chamber as the background value to subtract.
correctionMode = 'median';

% turn off the warning the Matlab modifies headers that contain
% unacceptable characters (like '.') to decluter the command line.
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%% loop through subfolders and load nuclei data into one table nr
nr = [];
for i=1:numel(dFolders)
%for i=1:1
    curFolder = dFolders(i).name;
    csvFolder = fullfile(fullfile(folderName,curFolder),'csvFIJI');
    analysisFileFound = isFileFound(csvFolder,'nucGeom.csv');

    if analysisFileFound == 0
        disp(['Could not find analysis file in subfolder ',csvFolder,'; excluding from nuclei analysis.']);
    else
        disp(['processing nuclei data from subfolder ',csvFolder,' ...']);
        
        % load nuclei geometry table 
        % coordinates in this file are relative to the *global image*
        ng = loadGlobalNucleiGeometry(csvFolder);
        
        % insert columns into nuclei geometry table that hold
        % - the eggChamber index, 
        % - the name of the subfolder holding the current eggChamber data 
        % - the name of the overarching folder holding all the eggChambers being analyzed
        ng = insertFolderInfoIntoTable(ng,i,folderName,curFolder);

        % compute # of nuclei
        nNuclei = size(ng,1); 
        disp(['Found ',num2str(nNuclei),' nuclei in subfolder ',curFolder,'...']);
        
        % load local coordinates, i.e. relative to the cropped nuclei
        % images. These will be used to transform cluster coordinates from
        % the local reference (cropped image) to a global reference (egg
        % chamber).
        ng = loadLocalNucleiGeometry(ng,nNuclei,csvFolder);

        curNr = ng;
        if nNuclei>0
            % compute # of channels
            [nc,allFilesPresent] = ...
                computeNumberOfChannelsInSubfolder(csvFolder);

            if allFilesPresent
                % loop through channels 
                for j=1:nc
                    % load nuclei intensity data channel
                    nucIntFileName = fullfile(csvFolder,...
                        ['C',num2str(j),'_nucInt.csv']);                 
                    nucIntRaw = readtable(nucIntFileName);

                    % append raw nuclei intensity data to nuclei table
                    curNr = addRawIntToTable(curNr,nucIntRaw,j,'nuc','_raw');

                    % load background intensity data channel 1 (whole-image) 
                    wholeImgIntFileName = fullfile(csvFolder,...
                        ['C',num2str(j),'_wholeImgInt.csv']);
                    wholeImgInt = readtable(wholeImgIntFileName);

                    % load background intensity data for the channel (egg chamber)
                    eggChamberIntFileName = fullfile(csvFolder,...
                        ['C',num2str(j),'_eggChamberInt.csv']);
                    eggChamberInt = readtable(eggChamberIntFileName);

                    % compute and append background-corrected nuclei intensity data (whole-image background)
                    if j==1
                        insertPosition = 'nucVolume';
                    else
                        insertPosition = ['nucC',num2str(j-1),'_Mode_wholeImgCorr'];
                    end
                    curNr = addImgCorrIntToTable(curNr,nucIntRaw,wholeImgInt,j,...
                        correctionMode,'nuc','wholeImgCorr',insertPosition,1);

                    % compute and append background-corrected nuclei intensity data (egg chamber background)
                    if j==1
                        insertPosition = 'nucVolume';
                    else
                        insertPosition = ['nucC',num2str(j-1),'_Mode_eggChamberCorr'];
                    end
                    curNr = addImgCorrIntToTable(curNr,nucIntRaw,eggChamberInt,j,...
                        correctionMode,'nuc','eggChamberCorr',insertPosition,1);

                    % load nucleoplasms geometry metrics (Volume)
                    if j == 1
                        plasmGeom = loadPlasmGeom(csvFolder,nNuclei);
                        for k=1:numel(plasmGeom.Properties.VariableNames)
                            oldVarName = plasmGeom.Properties.VariableNames{k};
                            newVarName = ['plasm',oldVarName];
                            curNr = addvars(curNr,plasmGeom.(oldVarName),...
                                'NewVariableNames', newVarName,...
                                'Before','rootFolder');
                        end
                    end

                    % load nucleoplasms intensity data channel
                    plasmIntRaw = loadPlasmInt(csvFolder,j,nNuclei);

                    % append raw nucleoplasm intensity data to nuclei table
                    curNr = addRawIntToTable(curNr,plasmIntRaw,j,'plasm','_raw');
                    
                    % compute and append background-corrected nucleoplasm intensity data (whole-image background)
                    if j==1
                        insertPosition = 'plasmVolume';
                    else
                        insertPosition = ['plasmC',num2str(j-1),'_Mode_wholeImgCorr'];
                    end
                    curNr = addImgCorrIntToTable(curNr,plasmIntRaw,wholeImgInt,j,...
                        correctionMode,'plasm','wholeImgCorr',insertPosition,0);

                    % compute and append background-corrected nucleoplasm intensity data (egg chamber background)
                    if j==1
                        insertPosition = 'plasmVolume';
                    else
                        insertPosition = ['plasmC',num2str(j-1),'_Mode_eggChamberCorr'];
                    end
                    curNr = addImgCorrIntToTable(curNr,plasmIntRaw,eggChamberInt,j,...
                        correctionMode,'plasm','eggChamberCorr',insertPosition,0);
                end
            end
        end
        if isempty(nr)
            nr = curNr;
        else
            % append current eggchamber data to full dataset table
            nr = [nr;curNr];
        end
    end
end

% add unique nucleus label
nr = addvars(nr,(1:size(nr,1))','NewVariableNames',{'nucUniqueID'},'Before','eggChamberID');

%% loop through subfolders and load cluster data
cl = [];
for i=1:numel(dFolders)
%for i=1:1
    curFolder = dFolders(i).name;
    csvFolder = fullfile(fullfile(folderName,curFolder),'csvFIJI');
    analysisFileFound = isFileFound(csvFolder,'nucGeom.csv');

    if analysisFileFound == 0
        disp(['Could not find analysis file in subfolder ',csvFolder,'; excluding from cluster analysis.']);
    else
        disp(['processing cluster data subfolder ',csvFolder,'...']);

        % compute # of nuclei
        ng = loadGlobalNucleiGeometry(csvFolder);
        nNuclei = size(ng,1); 
        disp(['Found ',num2str(nNuclei),' nuclei in subfolder ',csvFolder,'...']);

        % load local coordinates, i.e. relative to the cropped nuclei
        % images. These will be used to transform cluster coordinates from
        % the local reference (cropped image) to a global reference (egg
        % chamber).
        ng = loadLocalNucleiGeometry(ng,nNuclei,csvFolder);

        % loop through nuclei
        if nNuclei>0

            % compute # of channels
            [nc,allFilesPresent] = ...
                computeNumberOfChannelsInSubfolder(csvFolder);

            if allFilesPresent

                for k = 1:nNuclei
                
                    % load cluster geometry from Channel 1, plus folder
                    % name info
                    uniqueID = getNucleusUniqueID(nr,i,k);
                    curCl = loadClusterGeom(folderName,curFolder,i,k,uniqueID);

                    if ~isempty(curCl)
    
                        % loop through channels 
                        for j=1:nc
                        
                            % load raw cluster intensity data
                            fileName = ['C',num2str(j),'_nuc',num2str(k),'_clustInt_raw.csv'];
                            rawInt = readtable(fullfile(csvFolder,fileName));
    
                            % append raw cluster intensity data to cluster table
                            curCl = addRawIntToTable(curCl,rawInt,j,'clust','_raw');

                            % load nucleoplasm-corrected cluster intensity data
                            fileName = ['C',num2str(j),'_nuc',num2str(k),'_clustInt_plasmCorr.csv'];
                            plasmInt = readtable(fullfile(csvFolder,fileName));
    
                            % append nucleoplasm-corrected cluster intensity data to cluster table
                            curCl = addRawIntToTable(curCl,plasmInt,j,'clust','_plasmCorr');
    
                            % load nucleoplasm background intensity data for the channel
                            fileName = ['C',num2str(j),'_nuc',num2str(k),'_plasmInt.csv'];
                            plasmInt = readtable(fullfile(csvFolder,fileName));
    
                            % append columns containing nucleoplasm intensity data
                            curCl = addNucleoplasmIntensityToClusterTable(curCl,plasmInt,j);
 
                        end

                        % using nuclei local coordinates (relative to cropped image around nucleus) 
                        % & global coordinates(relative to entire egg chamber), 
                        % compute cluster global coordinates.
                        curCl = addClusterGlobalCoordinatesToTable(curCl,ng,k);
            
                        % append current nucleus clusters data to global table
                        if isempty(cl)
                            cl = curCl;
                        else
                            cl = [cl;curCl];
                        end
                    end
                end
            end
        end
    end
end

%% add cluster stats to nucleus table

% loop through egg chamber nuclei
eggChamberIDs = unique(nr.eggChamberID);
for i=1:numel(eggChamberIDs)
    idx = nr.eggChamberID == eggChamberIDs(i);
    idx1 = find(idx);
    idx1 = idx1(1);

    % compute # of channels
    csvFolder = fullfile( ...
        fullfile(nr.rootFolder(idx1),nr.eggChamberFolder(idx1)),'csvFIJI' );
    [nc,allFilesPresent] = computeNumberOfChannelsInSubfolder(csvFolder);
    
    % loop through nuclei and add summary of cluster metrics for each
    % nucleus to nr table
    nucIDs = unique(nr.nucID);
    for j=1:numel(nucIDs)
        nr = addSummaryClusterMetricsToNucleusTable(nr,cl,i,j,nc,correctionMode);
    end
end

%% add nuclei stats to cluster table

eggChamberIDs = unique(nr.eggChamberID);
for i=1:numel(eggChamberIDs)
    idx = nr.eggChamberID == eggChamberIDs(i);
    idx1 = find(idx);
    idx1 = idx1(1);

    % compute # of channels
    [nc,allFilesPresent] = computeNumberOfChannelsInSubfolder(...
        fullfile(fullfile(nr.rootFolder(idx1),nr.eggChamberFolder(idx1)),'csvFIJI'));
    
    % loop through nuclei and add summary of cluster metrics for each
    % nucleus to nr table
    nucIDs = unique(cl.nucID);
    for j=1:numel(nucIDs)
        idxN = (nr.eggChamberID == eggChamberIDs(i)) ...
            & (nr.nucID == nucIDs(j));

        if nr.nClusters(idxN) ~=0
            cl = addNucleiMetricsToClusterTable(...
                nr,cl,eggChamberIDs(i),nucIDs(j),nc);
        end
    end
end

%% save full data tables
writetable(nr,fullfile(folderName,'nucleusDataFull.txt'),'Delimiter','\t');
writetable(cl,fullfile(folderName,'clusterDataFull.txt'),'Delimiter','\t');

%% extract simplified nucleus table
% compute # of channels
[nc,allFilesPresent] = ...
    computeNumberOfChannelsInSubfolder(csvFolder);
nv = getNucleiVariablesToKeepInSummaryTable(nc);

nrs = extractVarsFromTable(nr,nv);

%% extract simplified cluster table
[nc,allFilesPresent] = ...
    computeNumberOfChannelsInSubfolder(csvFolder);
cv = getClusterVariablesToKeepInSummaryTable(nc);

cls = extractVarsFromTable(cl,cv);

%% save streamlined data tables
writetable(nrs,fullfile(folderName,'nucleusDataShort.txt'),'Delimiter','\t');
writetable(cls,fullfile(folderName,'clusterDataShort.txt'),'Delimiter','\t');