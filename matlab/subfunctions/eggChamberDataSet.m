
classdef eggChamberDataSet < handle

    properties (GetAccess = 'public', SetAccess = 'public')
        
        % table holding the data relative to each nucleus
        nucT = table();

        % exhaustive dataset before removing useless variables and metrics
        % left in case in-depth QC is needed.
        nucFullT = table();

        % table holding the  data relative to individual clusters
        clustT = table();
        
        % list of channels, e.g. [1,2,3,4]
        channelList = [];
        
        % number of Conditions
        nConditions = []; 

        % list of condition indices (num array) and corresponding folder names (cell array) 
        % condIndices is a nConditions x 1 numerical array which holds the
            % index of each condition  (usually will be 1:nConditions in a vertical vector)
            % e.g. [1;2]
        condIndices = [];
        % conditionNames is a nConditions x 1 cell array where the ith
            % entry is the name of the ith condition
            % e.g. {'Ctrl';
            %       'TRI';}
        conditionNames = {};

        % list of sample indices (num array) and the number of samples per condition
        % nSamples is an nConditions x 1 numerical array where ith entry 
            % is the number of samples in condition i).
            % e.g. [3;3]
        nSamples = [];
        % sampleIndices is an nConditions x 1 cell array where ith entry is
            % an nSamples(i) x 1 numeric array holding the index of each
            % sample (usually will be 1:nSamples(i) in a vertical vector)
            % e.g. sampleIndices{2,1} = [1;2;3]
        sampleIndices = {};
        % sampleNames is an nConditions x 1 cell array where ith entry is
            % a cell array (nSamples(i) x1) where the jth entry is the name
            % of the sample
            % e.g. sampleNames{2,1} = {'9-TRI-646,MPM2-488,Ser5ph-Cy3-1001zCorr';
            %                          '9-TRI-646,MPM2-488,Ser5ph-Cy3-1002zCorr';
            %                          '9-TRI-646,MPM2-488,Ser5ph-Cy3-1zCorr';}
        sampleNames = {}; 

        % input Folder path, e.g. '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img'
        inFolder = '';

        % index of the data channel that holds the egg chamber ID (scalar; set to zero if no
            % such info available).
            % e.g. 5
        eggChamberSegChannel = 0;
        
        % list of eggchamber IDs and their corresponding stages
        % eggChamberNumber is a nested cell array where eggChamberNumber{i}{j} is the number of egg
            % chambers segmented in condition i, sample j.
        eggChamberNumber = {};
        % eggChamberIDs is a nested cell array where eggChamberIDs{i}{j} is a eggChamberNumber{i}{j} x 1
            % numerical array holding the IDs of the segmented egg chambers in condition i, sample j.
        eggChamberIDs = {};
        % eggChamberStages is a nested cell array where eggChamberStages{i}{j} is a eggChamberNumber{i}{j} x 1
            % numerical array holding the stages of the segmented egg chambers in condition i, sample j.
            % defaulted to zero if the egg chamber stage info is missing.
        eggChamberStages = {};
        % eggChamberNumNucPerEC is a nested cell array where eggChamberNumNucPerEC{i}{j} is a eggChamberNumber{i}{j} x 1
            % numerical array holding the number of nuclei per segmented egg chambers in condition i sample j.
        eggChamberNumNucPerEC = {};
        
        % minimum volume of clusters in um^3 to be considered bona fide
        % HLB.
        clusterMinVol = 1;
    end

    properties (GetAccess = 'private', SetAccess = 'private')

        % ******************** eggchamber data I/O properties

        % name of the subfolder where the FOV-level (nucleus-wise, sampleROI-wise and wholeImg-wise) data is stored
        sampleFOVCsvFolderName = 'eggChamberCSV';
        % The folder should contain following files:
        % - allNucGeom.csv: geometry metrics of each nucleus
        % - C1_allNucInt.csv: nuclear-wise intensity metrics in channel 1
        % - C1_sampleROIInt.csv: sampleROI-wise intensity metrics in channel 1 (sample ROI is one giant ROI that covers broadly all egg chambers of the image)
        % - C1_wholeImgInt.csv: image-wise intensity metrics in channel 1
%         sampleFOVWiseFileNames = {'allNucGeom.csv';...
%             '_allNucInt.csv';...
%             '_wholeImgInt.csv';...
%             '_sampleROIInt.csv'};
        sampleFOVWiseFileNames = {'allNucGeom.csv';...
            '_allNucInt.csv';...
            '_wholeImgInt.csv';...
            '_sampleROIInt.csv';...
            'allEggChambersGeom.csv';...
            '_allEggChambersInt.csv';};

        % whether each file name of the previous list is channel-dependent or not
        sampleFOVWiseChannelDependent = logical([0;...
            1;...
            1;...
            1;...
            0;...
            1;]);
    
        % whether each file name of the previous list contains a single row (whole image data) or
        % one row per nucleus
        % 0 = single row table; 1 = 1 row per nucleus; 2 = 1 row per egg
        % chamber
        sampleFOVWiseSingleRow = [0;...
            0;...
            1;...
            1;...
            2;...
            2; ];
        
        % prefix to add in front of each variable name when loading the
        % corresponding table
        sampleFOVWisePrefix ={'nuc';...
            'nuc';...
            'wholeImg';...
            'sampleROI';...
            'eggChamber';
            'eggChamber';};

        % name of the optional subfolder where the 2D segmentations of the egg chambers are stored
        eggChamberSegFolderName = 'eggChamberSEG';
        % if used, the folder should hold a file called:
        eggChamberSegFileName = 'eggChamberStages.csv';
    
        % suffix to add at the end of each variable name when loading the
        % corresponding table
        eggChamberWiseSuffix ={'';...
            'raw';...
            'raw';...
            'raw';...
            '';...
            'raw';};
        
        % maximum number of channels (used when searching for data files) 
        maxChannels = 10;

        % ******************** nucleus & cluster I/O properties
        %  name of the subfolder where the cluster data is stored
        clusterFolderName = 'nucCSV';
        
        % file names (of name base in the case of channel dependent file names)
        clusterWiseFileNames = {'_plasmGeom.csv';...
            '_nucleoliGeom.csv';...
            '_plasmInt.csv';...
            '_nucleoliInt.csv';...
            '_clustInt_raw.csv';...
            '_clustInt_plasmCorr.csv'};

        % whether each file name of the previous list is channel dependent or not
        clusterWiseChannelDependent = logical([0;0;1;1;1;1;]);

        % whether each file name of the previous list contains a single row (whole image data) or
        % one row per nucleus
        clusterWiseSingleRow = logical([1;1;1;1;0;0;]);

        % prefix to add in front of each variable name when loading the
        % corresponding table
        clusterWisePrefix ={'plasm';'nucleoli';'plasm';'nucleoli';'clust';'clust';};

        % suffix to add at the end of each variable name when loading the
        % corresponding table
        clusterWiseSuffix ={'';'';'raw';'raw';'raw';'plasmCorr';};

        % ******************** Nuclei Table variables properties

        % the variables (i.e. columns) which store the sample and condition ID
        sampleVarList = {'cond_Idx', 'sample_Idx','sample_InputFileName'};
        
        % prefixes 
        % in a variable name, the prefix indicates what structure/ROI type the metric is applied to
        nucPrefix = 'nuc';
        wholeImgPrefix = 'wholeImg';
        sampleROIPrefix = 'sampleROI';
        nucleoliPrefix = 'nucleoli';
        plasmPrefix = 'plasm';
        clusterPrefix = 'clust';
        eggChamberPrefix = 'eggChamber';
        % the following 4 are obsolete as cluster summary metrics now
        % include the min vol and max vol values in their names
        % corresponding prefixes are added dynamically when the summary
        % metrics are calculated in addAverageClusterStatsToNucTable
        nucAvgClustPrefix = 'nucAvgClust'; % 
        nucStdClustPrefix = 'nucStdClust'; % 
        nucAvgClustMinVolPrefix = 'nucAvgClustMinVol'; 
        nucStdClustMinVolPrefix = 'nucStdClustMinVol'; 
        prefixList;
    
        % list of the variable basenames that pertain to each nucleus and its geometry - note
        % that the eggChamber_Idx and eggChamber_Stage fall into that category
        geomVarsBaseNameList = {'Label', 'Volume', 'SurfaceArea','Idx','Stage',......
            'MeanBreadth', 'Sphericity', 'EulerNumber', ...
        'BoxXMin', 'BoxXMax', 'BoxYMin', 'BoxYMax', ...
        'BoxZMin', 'BoxZMax', 'CentroidX', 'CentroidY', 'CentroidZ', ...
        'Dist',...
        'ElliCenterX','ElliCenterY','ElliCenterZ',...
        'InscrBallCenterX','InscrBallCenterY','InscrBallCenterZ','InscrBallRadius'...
        'ElliR1','ElliR2','ElliR3','ElliR1R2','ElliR2R3','ElliR1R3',...
        'ElliAzim','ElliElev','ElliRoll','VoxelCount',...
        'NumClusters','NumClustersMinVol'}; 
        % note 'NumClusters','NumClustersMinVol' are obsolete as cluster summary metrics now
        % include the min vol and max vol values in their names
        % corresponding baseNames are added dynamically when the summary
        % metrics are calculated in addAverageClusterStatsToNucTable
        % make sure that no entry in geomVarsBaseNameList is also part of the name of one of the entries in channelVarsBaseNameList
    
        % list of the variable basenames that pertain to an intensity metric
        channelVarsBaseNameList = {'Mean', 'StdDev', 'Max', 'Min','Median',...
            'Mode', 'Skewness', 'Kurtosis','Volume','NumberOfVoxels'}; 
        % the fact that 'Mean' is also part of a variable ('MeanBreadth') in geomVarsBaseNameList shouldnt be an issue
        % note that Volume is both an intensity and a geometry metric,
        % annoyingly

        % name of variables to remove in streamlined table
        geomVarsToRemoveInStreamLinedTable = {'SurfaceArea','MeanBreadth', ...
        'Sphericity', 'EulerNumber', ...
        'BoxXMin', 'BoxXMax', 'BoxYMin', 'BoxYMax', ...
        'BoxZMin', 'BoxZMax', 'CentroidX', 'CentroidY', 'CentroidZ', ...
        'Dist',...
        'ElliCenterX','ElliCenterY','ElliCenterZ',...
        'InscrBallCenterX','InscrBallCenterY','InscrBallCenterZ','InscrBallRadius'...
        'ElliR1','ElliR2','ElliR3','ElliR1R2','ElliR2R3','ElliR1R3',...
        'ElliAzim','ElliElev','ElliRoll',...
        'VoxelCount'};
        channelVarsToRemoveInStreamLinedTable = {'Min','Max','Mode', ...
            'Skewness', 'Kurtosis','Volume','NumberOfVoxels'};
    
        % suffixes
        % in an intensity variable name, the suffix indicates the processing applied to
        % the intensity
        rawSuffix = 'raw';
        plasmCorrSuffix = 'plasmCorr';
        eggChamberSubtractedSuffix = 'eggChamberSubtracted';
        sampleROISubtractedSuffix = 'sampleROISubtracted';
        wholeImgSubtractedSuffix = 'wholeImgSubtracted';
        nucleoliSubtractedPrefix = 'nucleoliSubtracted';
        plasmSubtractedPrefix = 'plasmSubtracted';
        suffixList;

    
        % variable names are built as follows:
            % sample ID variables: the full name is in the list sampleVarList
    
            % geometry variables: <prefix><geomVarBaseName>
    
            % channel variables:
            % <prefix>_C<channelIdx><channelVarsBaseNameList>_<suffix>
        
        % nucleus intensity background subtraction settings
        % these variables will get the bg intensity subtracted from them
        varsToBeSubtracted = {...
                'Mean',...
                'Max',...
                'Min',...
                'Median',...
                'Mode',...
            };
        
        % lits of the prefixes of the ROIs to use as reference for background subtraction
        % of nuclei intensities (wholeImg and eggChamber)
        nucBackgroundIntensityPrefixList;

        % lits of the prefixes of the ROIs to use as reference for background subtraction
        % of cluster intensities (plasm and nucleoli)
        clustBackgroundIntensityPrefixList;

        % when performing background subtraction, use the median intensity
        % of the reference region.
        metricToSubtract = 'Median';

        % nucleus variables to include in cluster table
        nucVarsToIncludeInClustTable = {'sample_InputFileName','cond_Idx',...
                'sample_Idx','eggChamber_Idx','eggChamber_Stage','nuc_Label'};

        % ************************** plotting settings    

        % used to set the spacing between cloud of points belonging to different nuclei    
        spacingUnit;

        % how much separation there is between conditions in units of spacingUnit
        condSeparator; 

        % how much separation there is between samples in units of spacingUnit
        sampleSeparator;

        % how much separation there is between egg chambers in units of spacingUnit
        ecSeparator;

        % width of the cloud of points for each egg chamber
        cloudWidth;
        
    end

    methods (Access = 'public')

        %% initialize object from input root folder
        function obj = eggChamberDataSet(inputFolder)
            obj.inFolder = inputFolder;
            obj.collectConditionsAndSamples;
            obj.resetPlottingDistances;

            obj.prefixList = {obj.nucPrefix,obj.wholeImgPrefix,...
                obj.sampleROIPrefix,obj.clusterPrefix,obj.nucleoliPrefix,...
                obj.plasmPrefix,obj.eggChamberPrefix,...
                obj.nucAvgClustMinVolPrefix,obj.nucStdClustMinVolPrefix;};
            obj.suffixList = {obj.rawSuffix,obj.plasmCorrSuffix,obj.eggChamberSubtractedSuffix,...
                obj.sampleROISubtractedSuffix,obj.wholeImgSubtractedSuffix,...
                obj.nucleoliSubtractedPrefix,obj.plasmSubtractedPrefix};
            obj.nucBackgroundIntensityPrefixList = {obj.wholeImgPrefix,obj.sampleROIPrefix,obj.eggChamberPrefix,obj.nucleoliPrefix};
            obj.clustBackgroundIntensityPrefixList = {obj.plasmPrefix,obj.nucleoliPrefix};

        end

        %% resets the plotting distances
        function resetPlottingDistances(obj)
            % used to set the spacing between cloud of points belonging to different nuclei    
            obj.spacingUnit = 1;
    
            % how much separation there is between conditions in units of spacingUnit
            obj.condSeparator = 2; 
    
            % how much separation there is between samples in units of spacingUnit
            obj.sampleSeparator = 2;
    
            % how much separation there is between egg chambers in units of spacingUnit
            obj.ecSeparator = 1;
    
            % width of the cloud of points for each egg chamber
            obj.cloudWidth = 0.3;
        end

        %%
        function setPlottingDistances(obj,userDefinedSpacingUnit,userDefinedCondSeparator,...
                userDefinedSampleSeparator,userDefinedEcSeparator,userDefinedCloudWidth)
            % used to set the spacing between cloud of points belonging to different nuclei    
            obj.spacingUnit = userDefinedSpacingUnit;
    
            % how much separation there is between conditions in units of spacingUnit
            obj.condSeparator = userDefinedCondSeparator; 
    
            % how much separation there is between samples in units of spacingUnit
            obj.sampleSeparator = userDefinedSampleSeparator;
    
            % how much separation there is between egg chambers in units of spacingUnit
            obj.ecSeparator = userDefinedEcSeparator;
    
            % width of the cloud of points for each egg chamber
            obj.cloudWidth = userDefinedCloudWidth;
        end
        
        %% load  nuclei data from all conditions into nuc data table 
        % this will look for the egg chamber stage info if present
        % (this will not load cluster data).
        function loadAllEggChamberNucleiData(obj)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames'); % mute unneeded warnings
            removeNeighbors = 1; % exclude from the data import the variables with neighbor in their name - we do not use them.
            obj.nucT = [];
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    disp(['Loading data from condition ',obj.conditionNames{i},...
                        ', sample ',obj.sampleNames{i}{j},' ...']);
                    obj.nucT = obj.combineEcTables(obj.nucT, loadEggChamberData(obj,i,j,removeNeighbors));
                end
            end
            if obj.nConditions == 0 || sum(obj.nSamples) == 0
                disp('Dataset is empty, cannot load data. Verify paths and initialize eggChamber object again.');
                return
            end
            obj.nucFullT = obj.nucT; % nucFullT backup copy of the exhaustive imported data table
            obj.getChannelList;
            obj.getEggChamberIDs;
            disp('done.');
        end    

        %% load cluster wise data from all conditions into nuc and clust data tables
        % populates the table clustT with cluster metrics
        % also adds nucleoli and nucleoplasm metrics to nucT and nucFullT tables
        % does NOT add summary cluster metrics to nucT and nucFullT tables
        function loadAllClusterData(obj)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames'); % mute unneeded warnings
            removeNeighbors = 1; % exclude from the data import the variables with neighbor in their name - we do not use them.
            
            % loop through conditions, find nuclei and load corresponding
            % data
            datasetNucTable = table();
            datasetNucFullTable = table();
            datasetClustTable = table();
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    
                    % find number/IDs of nuclei expected in the current sample based
                    % on nuclei table
                    nucIDList = obj.getExpectedNucIDsFromTable(i,j);
                    
                    % loop through nuclei and load the data
                    sampleNucTable = table();
                    sampleClustTable = table();
                    sampleNucFullTable = table();

                    for k=1:numel(nucIDList)
                        disp(['Loading data from condition ',obj.conditionNames{i},...
                        ', sample ',obj.sampleNames{i}{j},...
                        ', nucleus ',num2str(nucIDList(k)),' ...']);

                        % generate nuc & cluster table for current nucleus 
                        % that contain all the data that is available.
                        % curNucT table is nucleus metrics (nucleoplasm and nucleoli metrics)
                        % curClustT is cluster metrics.
                        [curNucT,curClustT] = obj.loadSampleClustTables(...
                            i,j,nucIDList(k),removeNeighbors);

                        % apend nucleus-wide info to the left of the cluster table
                        if ~isempty(curClustT)
                            curClustT = ...
                                obj.appendNucVarsToClustTable(curClustT,i,j,nucIDList(k)); 
                            
                        end

                        % append each of the newly created tables to current sample nuc and clust tables
                        sampleNucTable = obj.combineEcTables(sampleNucTable, curNucT);
                        sampleNucFullTable = obj.combineEcTables(sampleNucFullTable, curNucT);
                        sampleClustTable = obj.combineEcTables(sampleClustTable, curClustT);
                    end

                    % join table holding new nuc variables with existing nuc table for that sample
                    % using the nuclei labels as keys
                    sampleNucTable = join(obj.nucT(...
                        obj.nucT.cond_Idx == i & obj.nucT.sample_Idx == j,:),...
                        sampleNucTable,'Keys','nuc_Label');

                    sampleNucFullTable = join(obj.nucFullT(...
                        obj.nucFullT.cond_Idx == i & obj.nucFullT.sample_Idx == j,:),...
                        sampleNucFullTable,'Keys','nuc_Label');
                        
                    % append sample tables to dataset tables 
                    datasetNucTable = obj.combineEcTables(datasetNucTable,sampleNucTable);
                    datasetNucFullTable = obj.combineEcTables(datasetNucFullTable,sampleNucFullTable);
                    datasetClustTable = obj.combineEcTables(datasetClustTable,sampleClustTable);
                end
            end

            % update object tables nucT, fullT, clustT 
            obj.nucT = datasetNucTable;
            obj.nucFullT = datasetNucFullTable;
            obj.clustT = datasetClustTable;
            
            % streamline tables
            obj.streamLineTable('nuc');
            obj.streamLineTable('clust');
            disp('Done.');
        end

        %% add nucleus stats to cluster table
        % collects all variables from nucT table, excluding any that are
        % duplicates with cluster table variables, and adds then as new variables on the
        % cluster table Each row of the cluster table - i.e. each cluster - gets assigned 
        % values of the metric for the nucleus it belongs to.
        % (nuc_Label column is used to match clusters with their parent
        % nuclei).
        function tJoin = addNucStatsToClustTable(obj)
            % select only the variables that aren't already in clusters
            nucStatsVars = setdiff(obj.nucFullT.Properties.VariableNames, ...
                obj.clustT.Properties.VariableNames);
            nucStatsVars = [nucStatsVars,'nuc_Label'];
            idx = find(ismember(obj.nucFullT.Properties.VariableNames,nucStatsVars));
            
            n = obj.nucFullT(:,idx);

            % list of nuclei indices
            nucIdxRows = [obj.nucFullT.cond_Idx(:),obj.nucFullT.sample_Idx(:),obj.nucFullT.nuc_Label(:)];
            nucIdx = unique(nucIdxRows,'rows'); % nuclei shouldnt be repeated but whatever

            % list of indices in the cluster table
            clustIdxRows = [obj.clustT.cond_Idx(:),obj.clustT.sample_Idx(:),obj.clustT.nuc_Label(:)];
            
            % add columns in clustT table for nuc variables, all
            % initialized to either NaN or ''
            varsToAdd = setdiff(n.Properties.VariableNames,obj.clustT.Properties.VariableNames);
            idxVarsToAdd = ismember(n.Properties.VariableNames,varsToAdd);
            n2 = n(:,idxVarsToAdd);
            numericVars = varfun(@isnumeric,n2,'output','uniform');
            tJoin = obj.clustT;
            for i=1:size(n2,2)
                if numericVars(i)
                    tJoin = addvars(tJoin,NaN(size(tJoin,1),1),'NewVariableNames',n2.Properties.VariableNames(i));
                else
                    tJoin = addvars(tJoin,repmat({''},size(tJoin,1),1),'NewVariableNames',n2.Properties.VariableNames(i));
                end 
            end
            
            % loop through nuclei
            for ctr=1:size(nucIdx,1)
                % extract nuc stats for current nucleus
                nIdx = ismember(nucIdxRows,nucIdx(ctr,:),'rows');
                nT = n2(nIdx,:);
                nT = nT(1,:); % nuclei shouldnt be repeated but whatever
 
                % find rows matching the current nucleus in the cluster table
                cIdx = ismember(clustIdxRows,nucIdx(ctr,:),'rows');

                % add nuclei values to current rows
                tJoin(cIdx,end-size(n2,2)+1:end) = repmat(nT,sum(cIdx),1);
            end

            % make sure the variables are ordered
            cVarList = {'sample_InputFileName','cond_Idx','sample_Idx',...
                'eggChamber_Idx','eggChamber_Stage','nuc_Label','clust_Label','clust_Volume'};
            varOrder = 1:numel(cVarList);
            for i=1:numel(cVarList)
                if ismember(cVarList{i},tJoin.Properties.VariableNames)
                    tJoin = movevars(tJoin,cVarList{i},'Before',varOrder(i));
                end
            end
        end
        
        % remove nuc stats from cluster table
        function removeNucStatsToClustTable(obj)
            varsToKeep = setdiff(obj.clustT.Properties.VariableNames,obj.nucFullT.Properties.VariableNames);
            varsToKeep = [varsToKeep,'sample_InputFileName','cond_Idx','sample_Idx','eggChamber_Idx','eggChamber_Stage','nuc_Label'];      

            obj.clustT = obj.clustT(:,ismember(obj.clustT.Properties.VariableNames,varsToKeep));

        end
        %%
        function addAverageClusterStatsToNucTable(obj,minClustVolume,maxClustVolume)
            
            % variable name for the number of clusters
            numClustersVarName = ['nuc_NumClusters',num2str(minClustVolume),'_',num2str(maxClustVolume)];
            numClustersVarName = strrep(numClustersVarName,'.','pt');
            
            c = obj.clustT;
            n = obj.nucFullT; 
            
            % remove all the cluster values outside of the volume range
            idxClusters = c.clust_Volume >= minClustVolume...
                          & c.clust_Volume <= maxClustVolume;
            
            c = c(idxClusters,:);
            
            %% collect all the variables from the cluster table except those that were
            % lifted from the nucleus table, using prefixes as filters for varibale
            % names
            varsToExclude = ismember(c.Properties.VariableNames,{'clust_Label'});
            pList = {'nuc_','sampleROI_','wholeImg_','eggChamber_','eggChamber_','plasm_','nucleoli_'};
            for i=1:numel(pList)
                varsToExclude = varsToExclude  | ...
                    cellfun(@startsWith,...
                    c.Properties.VariableNames,repmat(pList(i),...
                    size(c.Properties.VariableNames)));
            end
            
            % retain the few variables that are important but have forbiden prefixes
            varsToExclude(ismember(c.Properties.VariableNames,{'nuc_Label'})) = 0;
            varsToExclude(ismember(c.Properties.VariableNames,{'eggChamber_Idx'})) = 0;
            varsToExclude(ismember(c.Properties.VariableNames,{'eggChamber_Stage'})) = 0;
            
            % filter out excluded variables from cluster table
            c = c(:,~varsToExclude);
            
            %% average each cluster metric across each nucleus
            c = grpstats(c,...
                            ["cond_Idx","sample_Idx","sample_InputFileName",...
                            "eggChamber_Idx","eggChamber_Stage","nuc_Label"],...
                            ["mean","std"]);
            
            % rename the variable holding the number of clusters in each nucleus "nClusters"
            c = renamevars(c,{'GroupCount'},numClustersVarName);
            
            % remove from the averaged table the new row names created by the grpstats function
            c.Properties.RowNames = {}; 
            
            %% replace prefixes mean_ and std_ by nucAvgClust<MinVol>_<MaxVol>_ and nucStdClust<MinVol>_<maxVol>_
            prefixList1 ={'mean_','std_'}; % prefix to replace in t (generated by grpstats)
            prefixList2 ={'Avg','Std'}; % corresponding name of the metric in new variable
            newPrefix = {};
            for i=1:numel(prefixList1)
                vIn = c.Properties.VariableNames(cellfun(@startsWith,c.Properties.VariableNames,repmat(prefixList1(i),...
                        size(c.Properties.VariableNames))));
                newPrefix{i} = strrep(['nuc',prefixList2{i},'Clust',num2str(minClustVolume),'_',num2str(maxClustVolume),'_'],'.','pt');
                vOut = cellfun(@strrep,vIn,repmat(prefixList1(i),size(vIn)),...
                    repmat(newPrefix(i),size(vIn)),...
                    'UniformOutput',0);
                c = renamevars(c, vIn,vOut);
            end
            
            %% add rows for nuclei which do not contain clusters
             
            % find unique combos of nuclei IDs present in nucleus table but absent in
            % cluster table
            % "cond_Idx","sample_Idx","sample_InputFileName","eggChamber_Idx","eggChamber_Stage","nuc_Label"
            idCols = ["cond_Idx","sample_Idx","sample_InputFileName","eggChamber_Idx","eggChamber_Stage","nuc_Label"];
            idxNucCols = ismember(n.Properties.VariableNames,idCols);
            idxClustCols = ismember(c.Properties.VariableNames,idCols);
            missingNucs = setdiff( unique( n(:,idxNucCols) ,'rows'),...
                unique( c(:,idxClustCols) ,'rows'));
            
            % generate a new table with the same column as the cluster table
            if isempty(c)
                newNucs = array2table(zeros(size(missingNucs,1),size(c,2)), 'VariableNames',c.Properties.VariableNames);
            else
                newNucs = repmat(c(1,:),size(missingNucs,1),1);
            end
            
            % initialize all values to NaN or '', except number of clusters set to 0
            numericVars = varfun(@isnumeric,newNucs,'output','uniform');
            newNucs{:,numericVars} = NaN;
            newNucs(:,~numericVars) = {''};
            newNucs.(numClustersVarName)(:) = 0;
            
            % enter the labels corresponding to the nuclei without clusters
            for i=1:numel(idCols)
                newNucs.(idCols{i}) = missingNucs.(idCols{i});
            end
            
            % combine the missing nuclei with t
            t2 = [c;newNucs];
            
            %% join the averaged cluster metrics with the nucleus metrics table
            % join the averaged metrics with the nucleus table
            newNucTable = join(n,t2,"Keys",...
                ["cond_Idx","sample_Idx","sample_InputFileName",...
                "eggChamber_Idx","eggChamber_Stage","nuc_Label"]);
            
            %% re-order key variables to the left for easier browsing of the table
            
            
            % make sure the variables are ordered
            cVarList = {'sample_InputFileName','cond_Idx','sample_Idx',...
                'eggChamber_Idx','eggChamber_Stage','nuc_Label','nuc_Volume'};
            varOrder = 1:numel(cVarList);
            for i=1:numel(cVarList)
                if ismember(cVarList{i},newNucTable.Properties.VariableNames)
                    newNucTable = movevars(newNucTable,cVarList{i},'Before',varOrder(i));
                end
            end

            % place the number of clusters in the nucleus right after the nuc_Volume:
            if ismember(numClustersVarName,newNucTable.Properties.VariableNames) ...
                    && ismember('nuc_Volume',newNucTable.Properties.VariableNames)
                newNucTable = movevars(newNucTable, {numClustersVarName},'After','nuc_Volume');
            end
            
            obj.nucFullT = newNucTable;

            % add new variable names to the list of acceptable basenames and prefixes
            obj.addClusterVariableNames(numClustersVarName,newPrefix)

        end

        % remove summary metrics of cluster stats from nuc table if needed
        function deleteAverageClusterStatsFromNucTable(obj,minClustVolume,maxClustVolume)
            % variable name for the number of clusters
            numClustersVarName = ['nuc_NumClusters',num2str(minClustVolume),'_',num2str(maxClustVolume)];
            numClustersVarName = strrep(numClustersVarName,'.','pt');
            
            if ismember(numClustersVarName,obj.nucFullT.Properties.VariableNames)
                newNuc = removevars(obj.nucFullT,numClustersVarName);
            else
                newNuc = obj.nucFullT;
                disp(['Variable ',numClustersVarName,' is absent already.']);
            end
            
            % finding other metrics
            cPrefixList ={'Avg','Std'}; % corresponding name of the metric in new variable
            newPrefix = {};
            idxToRemove = false(size(newNuc.Properties.VariableNames));
            for i=1:numel(cPrefixList)
                newPrefix{i} = strrep(['nuc',cPrefixList{i},'Clust',num2str(minClustVolume),'_',num2str(maxClustVolume),'_'],'.','pt'); 
                idxToRemove = idxToRemove | cellfun(@startsWith,newNuc.Properties.VariableNames,repmat(newPrefix(i),...
                        size(newNuc.Properties.VariableNames)));
            end
            newNuc = newNuc(:,~idxToRemove);
            obj.nucFullT = newNuc;
        end

        % add variable names (generated when computing cluster summary metrics by nucleus) 
        % to the list of acceptable basenames and prefixes
        function addClusterVariableNames(obj,numClustersVarName,newPrefix)

            % number of clusters as a geometric variable
            if startsWith(numClustersVarName,[obj.nucPrefix,'_'])
                numClustersVarName = numClustersVarName((strlength(obj.nucPrefix)+2):end);
            end
            obj.geomVarsBaseNameList = [obj.geomVarsBaseNameList,numClustersVarName];

            % cluster mean and Std as new prefixes
            for i=1:numel(newPrefix)
                obj.prefixList = [obj.prefixList,newPrefix{i}];
            end
        end

        %% add average cluster stats to nuc tables
        function addAverageClusterStatsToNucTableOld(obj)
            % generate list of cluster variables to average out
            cv = obj.clustT.Properties.VariableNames;
            idx = find(cell2mat( cellfun( @contains, cv,...
                        repmat({'clust'}, size(cv)),...
                        'UniformOutput',0) )); 
            cIdx = find(ismember(cv,{'clust_Label'}));
            idx = setdiff(idx,cIdx); % indices of the variables to average (in clustT table)            
            
            t = table();
            for ctr=1:size(obj.nucFullT,1)
                i = obj.nucFullT.cond_Idx(ctr);
                j = obj.nucFullT.sample_Idx(ctr);
                k = obj.nucFullT.nuc_Label(ctr);
                
                % extract cluster values for current nucleus (ALL clusters)
                cT1 = obj.clustT(...
                   obj.clustT.cond_Idx == i...
                   & obj.clustT.sample_Idx == j ...
                   & obj.clustT.nuc_Label == k,idx);

                cT1 = compileClusterStatsForCurrentNucleus(cT1,...
                    obj.nucAvgClustPrefix,obj.nucStdClustPrefix,...
                    [obj.nucPrefix,'_NumClusters']);
                
                % extract cluster values for current nucleus (ONLY clusters with volume above threshold) 
                cT2 = obj.clustT(...
                   obj.clustT.cond_Idx == i...
                   & obj.clustT.sample_Idx == j ...
                   & obj.clustT.nuc_Label == k...
                   &obj.clustT.clust_Volume >= obj.clusterMinVol,idx);

                cT2 = compileClusterStatsForCurrentNucleus(cT2,...
                    obj.nucAvgClustMinVolPrefix,obj.nucStdClustMinVolPrefix,...
                    [obj.nucPrefix,'_NumClustersMinVol']);

                cT = join(cT1,cT2,'Keys','nuc_Label');

                % append current nuclei cluters metrics at the bottom of global table
                t = obj.combineEcTables(t,cT);  
            end
            
            % join horizontally t (average cluster metrics per nucleus) with the nucT table
            % holding all other nuclei metrics
            if size(t,1) ~= size(obj.nucFullT,1)
                disp(['Cannot append average cluster values to nucleus table; ',...
                    'Nuc table has ',num2str(size(obj.nucFullT,1)),' rows while ',...
                    'Clust table has ',num2str(size(t,1)),' rows.']);
            else
                % add dummy column to each variable to use as 
                dummyKey = (1:size(obj.nucFullT,1))';
                nuc2 = addvars(obj.nucFullT,dummyKey,'NewVariableNames',{'Key'});
                t = addvars(t,dummyKey,'NewVariableNames',{'Key'});
                
                % add new cluster averaged values to nuc table.
                tJoin = join(nuc2,t,'Keys','Key');

                % remove key
                tJoin = removevars(tJoin,{'Key'});

                % QC check
                if sum(abs(tJoin.nuc_Label_nuc2 - tJoin.nuc_Label_t)) ~=0
                    failedIndices = find(abs(tJoin.nuc_Label_nuc2 - tJoin.nuc_Label_t)~=0)
                    disp(['QC failed, nuc_Labels of cluster metrics and nuc table do not match;',...
                        ' cannot add cluster metrics to nucT.']);
                else
                    
                    tJoin = removevars(tJoin,{'nuc_Label_t'});
                    tJoin = renamevars(tJoin,{'nuc_Label_nuc2'},{'nuc_Label'});
                end
            end
            obj.nucFullT = tJoin;


            %%
            function tOut = compileClusterStatsForCurrentNucleus(tIn,...
                    curNucAvgClustPrefix,curNucStdClustPrefix,curNumClusterPrefix)
                
                if ~isempty(tIn)
                    nC = size(tIn,1); % number of clusters in current nucleus

                    % average/std cluster metrics and keep only meaningful variables
                    cT1 = varfun(@mean,tIn,'InputVariables',@isnumeric);
                    cT2 = varfun(@std,tIn,'InputVariables',@isnumeric);
                    
                    % update variable names from mean_clust to avgClust
                    newVars = cellfun(@strrep,cT1.Properties.VariableNames,...
                        repmat({'mean_clust'},size(cT1.Properties.VariableNames)),...
                        repmat({curNucAvgClustPrefix},size(cT1.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT1 = renamevars(cT1,cT1.Properties.VariableNames,newVars);
                    
                    % update variable names from std_clust to stdClust
                    newVars = cellfun(@strrep,cT2.Properties.VariableNames,...
                        repmat({'std_clust'},size(cT2.Properties.VariableNames)),...
                        repmat({curNucStdClustPrefix},size(cT2.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT2 = renamevars(cT2,cT2.Properties.VariableNames,newVars);
                    
                    % add a nuc_Label variable to join tables (and
                    % validation later)
                    cT1 = addvars(cT1,k,'NewVariableNames',{'nuc_Label'});
                    cT2 = addvars(cT2,k,'NewVariableNames',{'nuc_Label'});
                    tOut = join(cT1,cT2,'Keys','nuc_Label');

                    % compute the number of clusters per nucleus 
                    tOut = addvars(tOut,nC,'NewVariableNames',{curNumClusterPrefix});   
                else
                    % no clusters found in current nucleus, generating a
                    % row of NaN values
                    tIn = obj.clustT(1,idx);
                    for cc = 1:size(tIn,2)
                        tIn{1,cc} = NaN;
                    end

                    % update variable names from clust to avgClust
                    newVars = cellfun(@strrep,tIn.Properties.VariableNames,...
                        repmat({'clust'},size(tIn.Properties.VariableNames)),...
                        repmat({curNucAvgClustPrefix},size(tIn.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT1 = renamevars(tIn,tIn.Properties.VariableNames,newVars);

                    % same for stdClust
                    newVars = cellfun(@strrep,tIn.Properties.VariableNames,...
                        repmat({'clust'},size(tIn.Properties.VariableNames)),...
                        repmat({curNucStdClustPrefix},size(tIn.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT2 = renamevars(tIn,tIn.Properties.VariableNames,newVars);

                    % add a nuc_Label variable (used for joinig and validation later)
                    cT1 = addvars(cT1,k,'NewVariableNames',{'nuc_Label'});
                    cT2 = addvars(cT2,k,'NewVariableNames',{'nuc_Label'});
                    tOut = join(cT1,cT2,'Keys','nuc_Label');

                    % compute the number of clusters per nucleus 
                    tOut = addvars(tOut,0,'NewVariableNames',{curNumClusterPrefix});

                end
            end
        end

        %% remove unlikely to be used variables from Nuclei table
        % tableType should be either 'nuc' or 'clust'
        function streamLineTable(obj,tableType)
            disp(['Streamlining ',tableType,' table...']);
            [c,nChannels] = obj.getChannelList;
            switch tableType
                case 'nuc'
                    t = obj.nucT;
                case 'clust'
                    t = obj.clustT;
            end

            if ismember('eggChamber_Idx',t.Properties.VariableNames) ...
                    && ismember('eggChamber_Stage',t.Properties.VariableNames) ...
                    && obj.eggChamberSegChannel ~= 0
                removeEggChamberSegChannelVars = 1;
            else
                removeEggChamberSegChannelVars = 0;
            end
            
            % compile list of variables to remove
            varsToRemoveList = {};
            for i=1:numel(obj.prefixList)
                for n=1:numel(obj.geomVarsToRemoveInStreamLinedTable)
                    newVar = buildVarName(obj,obj.prefixList{i},0,...
                        obj.geomVarsToRemoveInStreamLinedTable{n},'','geom');
                    varsToRemoveList = [varsToRemoveList; newVar];
                end

                for j=1:nChannels
                    for k=1:numel(obj.suffixList)
                        % if the current channel is the one holding the egg
                        % chamber segmentation channel metrics, remove all
                        % metrics (the egg chamber ID and stage are already saved in a different column.)
                        if removeEggChamberSegChannelVars ...
                                && c(j) == obj.eggChamberSegChannel
                            curVarList = obj.channelVarsBaseNameList;
                        else
                            curVarList = obj.channelVarsToRemoveInStreamLinedTable;
                        end
                        
                        for n=1:numel(curVarList)
                            newVar = obj.buildVarName(obj.prefixList{i},c(j),...
                                curVarList{n},obj.suffixList{k},'channel');
                            varsToRemoveList = [varsToRemoveList; newVar];
                        end
                    end
                end
            end

            % remove variable names that arent present
            varsToRemoveList = varsToRemoveList(...
                ismember(varsToRemoveList,t.Properties.VariableNames));

            t = removevars(t,varsToRemoveList);
            
            switch tableType
                case 'nuc'
                    obj.nucT = t;
                case 'clust'
                    obj.clustT = t;
            end

            % reorder variables
            obj.reOrderNucleiVariables();

            % sort rows
            obj.sortNucleiRowsByEggChamber();

            disp('Done.');
        end

        %% perform background subtraction on nuclei intensity values
         function backgroundCorrectNucIntensity(obj)
            
            obj.nucFullT = bgCorrTable(obj.nucFullT);

             function tOut = bgCorrTable(tIn)
                tOut = tIn;

                [c,nChannels] = obj.getChannelList;

                % find column indices of the variables that need subtracting
                varList = tIn.Properties.VariableNames;
                idx = [];
                for i=1:numel(obj.varsToBeSubtracted)
                    curIdx = find(cell2mat( cellfun( @contains, varList,...
                        repmat(obj.varsToBeSubtracted(i), size(varList)),...
                        'UniformOutput',0) ));
                    idx = [idx,curIdx];
                end
                
                % remove variables relative to background regions from the 'to be
                % subtracted' list
                for i=1:numel(obj.nucBackgroundIntensityPrefixList)
                    curIdx = find(cell2mat( cellfun( @contains, varList,...
                        repmat(obj.nucBackgroundIntensityPrefixList(i), size(varList)),...
                        'UniformOutput',0) ));
                    idx = setdiff(idx,curIdx);
                end

                varList = tIn.Properties.VariableNames(idx);

                % loop through variables to subtract (wholeImg and sampleROI)
                for v = 1:numel(obj.nucBackgroundIntensityPrefixList)
                    
                    % loop through color channels
                    for i=1:nChannels

                        % build variable to subtract
                        if strcmp(obj.nucBackgroundIntensityPrefixList{v},'sampleROI')
                                % use the mode of the intensity within the
                                % sample ROI as the background value
                                subVarName = obj.buildVarName(...
                                    obj.nucBackgroundIntensityPrefixList{v},c(i),...
                                    'Mode',obj.rawSuffix,'channel');
                                subVar = tIn.(subVarName);

                        elseif strcmp(obj.nucBackgroundIntensityPrefixList{v}, 'wholeImg')
                                % use the geometric mean between the min and the mode 
                                % of the intensity within the
                                % whole img as the background value (+1 added to avoid geometric means equal to zero)
                                min1 = obj.buildVarName(...
                                    obj.nucBackgroundIntensityPrefixList{v},c(i),...
                                    'Min',obj.rawSuffix,'channel');
                                mode1 = obj.buildVarName(...
                                    obj.nucBackgroundIntensityPrefixList{v},c(i),...
                                    'Mode',obj.rawSuffix,'channel');
                                subVar = geomean([tIn.(min1)+1,tIn.(mode1)],2);
                        else
                            % other variables, use the median value
                            subVarName = obj.buildVarName(obj.nucBackgroundIntensityPrefixList{v},c(i),...
                                'Median',obj.rawSuffix,'channel');
                            if ismember(subVarName,tIn.Properties.VariableNames)
                                subVar = tIn.(subVarName);
                            else
                               disp(['Warning: Could not find variable ',subVarName,' to background-subtract nuclei intensity. Skipping.']); 
                               subVar = [];
                            end
                        end
                        
                        if ~isempty (subVar)
                            % find variables in tIn that are relative to the current color channel
                            idxC = cell2mat( cellfun( @contains, varList,...
                                repmat({['_C',num2str(c(i))]}, size(varList)),...
                                'UniformOutput',0) );
                            
                            curVarList = varList(idxC);
                            for k =1:numel(curVarList)
                                newVarName = strrep(curVarList{k},...
                                    ['_',obj.rawSuffix],...
                                    ['_',obj.nucBackgroundIntensityPrefixList{v},'Subtracted']);
        
                                newVar = tIn.(curVarList{k}) - subVar;
    
                                % overwrite background corrected variable if it already exists,
                                % create a new var otherwise
                                if ismember(newVarName,tIn.Properties.VariableNames)
                                    tOut.(newVarName) = newVar;
                                else
                                    tOut = addvars(tOut,newVar,'NewVariableNames',newVarName);
                                end
                            end
                        end
                    end
                end
            end
         end

         %% perform background subtraction on cluster intensity values
         function backgroundCorrectClustIntensity(obj)
            
            tIn = obj.clustT;
            tOut = tIn;

            [c,nChannels] = obj.getChannelList;
            
            % find column indices of the variables that contain the clust_
            % prefix
            varList = tIn.Properties.VariableNames;
            idx = cell2mat( cellfun( @contains, varList,...
                    repmat({'clust_'}, size(varList)),...
                    'UniformOutput',0) );
            varList = varList(idx);

            % find column indices of the variables that need subtracting
            idx = [];
            for i=1:numel(obj.varsToBeSubtracted)
                curIdx = find(cell2mat( cellfun( @contains, varList,...
                    repmat(obj.varsToBeSubtracted(i), size(varList)),...
                    'UniformOutput',0) ));
                idx = [idx,curIdx];
            end
            varList = varList(idx);

            %ignore non-raw variables
            idx = cell2mat( cellfun( @contains, varList,...
                    repmat({'_raw'}, size(varList)),...
                    'UniformOutput',0) );
            varList = varList(idx);

            % loop through variables to subtract (plasm and nucleoli)
            for v = 1:numel(obj.clustBackgroundIntensityPrefixList)

                % loop through color channels
                for i=1:nChannels

                    % build variable to subtract
                    subVarName = obj.buildVarName(...
                                obj.clustBackgroundIntensityPrefixList{v},c(i),...
                                'Median',obj.rawSuffix,'channel');
                    subVar = tIn.(subVarName);

                    % find variables in tIn that are relative to the current color channel
                    idxC = cell2mat( cellfun( @contains, varList,...
                        repmat({['_C',num2str(c(i))]}, size(varList)),...
                        'UniformOutput',0) );
                    
                    curVarList = varList(idxC);
                    for k =1:numel(curVarList)
                        newVarName = strrep(curVarList{k},...
                            ['_',obj.rawSuffix],...
                            ['_',obj.clustBackgroundIntensityPrefixList{v},'Subtracted']);

                        newVar = tIn.(curVarList{k}) - subVar;

                        % overwrite background corrected variable if it already exists,
                        % create a new var otherwise
                        if ismember(newVarName,tIn.Properties.VariableNames)
                            tOut.(newVarName) = newVar;
                        else
                            tOut = addvars(tOut,newVar,'NewVariableNames',newVarName);
                        end
                    end
                end
            end
            obj.clustT = tOut;
        end

        %% generate table holding summary statistics for each egg chamber
        function sumT = generateEggChamberSummaryTable(obj)
            conditionName = {};
            sampleName = {};
            eggChamber_Idx = [];
            eggChamber_Stage = [];
            numberOfNucleiInEggChamber = [];
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        conditionName = ...
                            [conditionName; obj.conditionNames{i}];
                        sampleName = [sampleName; obj.sampleNames{i}{j}];
                        eggChamber_Idx = ...
                            [eggChamber_Idx; obj.eggChamberIDs{i}{j}(k)];
                        eggChamber_Stage = ...
                            [eggChamber_Stage; ...
                            obj.eggChamberStages{i}{j}(k)];
                        numberOfNucleiInEggChamber = ...
                            [numberOfNucleiInEggChamber; ...
                            obj.eggChamberNumNucPerEC{i}{j}(k)];
                    end
                end
            end
            sumT = table(conditionName,sampleName,eggChamber_Idx,...
                eggChamber_Stage,numberOfNucleiInEggChamber);
        end
        
        %% scatter plot a metric by sample
        function scatterPlotNucleiMetricBySample(obj,prefix,channel,baseName,suffix)

            varName = obj.buildVarName(prefix,channel,baseName,suffix,'geom');
            
            % build figure
            figure('Name',strrep(varName,'_',' '));
            hold;
            
            % collect the x values to plot each condition/sample at
            [xSampleVals,xSampleIDs] = obj.getSampleXValuesBySample();

            % collect the values of the metric for each condition/sample
            xPlot = [];
            yPlot = [];
            xSampleValsVec = [];
            xSampleIDsVec = {};
            for j=1:numel(obj.condIndices)
                for k=1:obj.nSamples(j)
                    s = obj.sampleIndices{j};

                    % collect values for the desired metric from all nuclei for the
                    % current sample/condition
                    x = obj.nucFullT.(varName)(...
                        obj.nucFullT.cond_Idx ==obj.condIndices(j) ...
                        & obj.nucFullT.sample_Idx == s(k));
            
                    % generate slightly offset x coordinates for each nucleus,
                    % centered around the sample X
                    nNuclei = size(x,1);
                    if nNuclei >1
                        % spacing between nuclei
                        nucSpacing = obj.spacingUnit ...
                            * obj.cloudWidth/(nNuclei-1);
            
                        % x coordinate for each nucleus of current condition/sample
                        curXPlot = xSampleVals{j}(k) - floor(nNuclei/2)*nucSpacing ...
                            + (0:(nNuclei-1))*nucSpacing;
            
                        % y coordinate for each nucleus of current condition/sample
                        curYPlot = obj.nucFullT.(varName)(...
                            obj.nucFullT.cond_Idx ==obj.condIndices(j) ...
                            & obj.nucFullT.sample_Idx == s(k))';
            
                    elseif nNuclei == 1
                        curXPlot = xSampleVals{j}(k);
                        curYPlot = obj.nucFullT.(varName)(...
                            obj.nucFullT.cond_Idx ==obj.condIndices(j) ...
                            & obj.nucFullT.sample_Idx == s(k))';
            
                    elseif nNuclei == 0
                        curXPlot = [];
                        curYPlot = [];
                    end
            
                    % append coordinates of current condition/sample to global list
                    xPlot = [xPlot,curXPlot];
                    yPlot = [yPlot,curYPlot];
            
                    xSampleValsVec = [xSampleValsVec,xSampleVals{j}(k)];
                    xSampleIDsVec = [xSampleIDsVec,xSampleIDs{j}{k}];
                end
            end
            plot(xPlot,yPlot,'o');
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylabel(strrep(baseName,'_','\_'));

        end

        %% scatter plot arbitrary metric by egg chamber
        % older version, use scatterPlotNucArbitraryMetricByEggChamber
        % instead
        function fh = scatterPlotArbitraryMetricByEggChamberOld(obj,figHandle,NucOrClustTable,yMetric,idx,eggChamberStagesToInclude)
            
            % if shorthand 'all' was used for cell stages to include, replace it by
            % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % figure out which table to use
            NucOrClustTable = lower(NucOrClustTable);
            switch NucOrClustTable
                case 'nuc'
                    t = obj.nucFullT;
                case 'clust'
                    t = obj.clustT;
            end

            % check that the metric and index entries have the right size that matches that of the table.
            fh = [];
            if numel(yMetric) ~= size(t,1) 
                disp(['Error: metric to plot does not have the size matching the ',NucOrClustTable,' table']);
                return
            end
            if numel(idx) ~= size(t,1)
                disp(['Error: index does not have the size matching the ',NucOrClustTable,' table']);
                return
            end

            % build figure, or assign fh handle to desired figure
            if strcmp(get(figHandle,'type'),'figure')
                fh = figHandle;
                set(0,'CurrentFigure',fh);
                if ~ishold
                    hold;
                end
            else
                fh = figure('Name',strrep(varName,'_',' '));
                hold;
            end    
            
            % collect the x values to plot each condition/sample at
            [xEggChamberVals, xEggChamberIDs, ~] = ...
                getSampleXValuesByEggChamberGatedByIdx(obj,eggChamberStagesToInclude,t,idx);

            % generate color maps
            nc = numel(obj.condIndices);
            cmData =cbrewer('qual', 'Set1', max(nc,3));
            
            % collect the values of the metric for each condition/sample
            
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            for i=1:numel(obj.condIndices)
                xPlot = [];
                yPlot = [];
                s = obj.sampleIndices{i};
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude)
                            % collect values for the desired metric from all nuclei for the
                            % current sample/condition
                            x = yMetric(...
                                t.cond_Idx ==obj.condIndices(i) ...
                                & t.sample_Idx == s(j) ...
                                & t.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k) ...
                                & idx);

                            % generate slightly offset x coordinates for each nucleus,
                            % centered around the sample X
                            nNuclei = size(x,1);
                            if nNuclei >1
                                % spacing between nuclei
                                nucSpacing = obj.spacingUnit ...
                                    * obj.cloudWidth / (nNuclei-1);
                    
                                % x coordinate for each nucleus of current condition/sample
                                curXPlot = xEggChamberVals{i}{j}(k) ...
                                    - floor(nNuclei/2)*nucSpacing ...
                                    + (0:(nNuclei-1))*nucSpacing;
                    
                                % y coordinate for each nucleus of current condition/sample
                                curYPlot = obj.nucFullT.(varName)(...
                                    obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucFullT.sample_Idx == s(j)...
                                    & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 1
                                curXPlot = xEggChamberVals{i}{j}(k);
                                curYPlot = obj.nucFullT.(varName)(...
                                    obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucFullT.sample_Idx == s(j)...
                                    & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 0
                                curXPlot = [];
                                curYPlot = [];
                            end
                    
                            % append coordinates of current condition/sample to global list
                            xPlot = [xPlot,curXPlot];
                            yPlot = [yPlot,curYPlot];

                            xErr = mean(curXPlot);
                            yErr = mean(curYPlot);
                            eErr = std(curYPlot);

                            xSampleValsVec = [xSampleValsVec,xEggChamberVals{i}{j}(k)];
                            xSampleIDsVec = [xSampleIDsVec,xEggChamberIDs{i}{j}{k}];

                            p = scatter(curXPlot,curYPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                            alpha(p,0.3);
                            errorbar(xErr,yErr,eErr,'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                        end
                    end
                end
                if ~isempty(yPlot)
                    yMax = max(yMax,max(yPlot(:)));
                    yMin = min(yMin,min(yPlot(:)));
                end
                
            end
            
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylim([min([0,1.2*yMin]),1.2*yMax]);
            ylabel(strrep(baseName,'_','\_'));
            grid on
        end    

        %% scatter plot a nucleus metric by egg chamber (old)
        % older version, use scatterPlotNucTableMetricByEggChamber
        % instead
        function scatterPlotNucleiMetricByEggChamberOld(obj,prefix,channel,baseName,suffix,eggChamberStagesToInclude)
            % prefix: any allowable prefix which marks the compartment the metric is calculated on e.g. 'nuc', or 'clust'
            % channel: intensity channel , e.g. 1. 
                % (Value is ignored if the metric is a geometry feature 
                % rather than an intensity feature.
            % baseName: any allowable metric basename, e.g. 'Mean' or 'Volume'
            % suffix: any allowable suffix which marks processing steps,
            % e.g. 'raw' or 'eggChamberCorr'
            % if shorthand 'all' was used for cell stages to include, replace it by
            % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % check that the metric is present in the data table.
            varName = obj.buildVarName(prefix,channel,baseName,suffix,'geom');
            if ~ismember( varName, obj.nucFullT.Properties.VariableNames)
                disp(['Variable ',varName,' absent from table, cannot plot.']);
                return
            end

            % build figure
            figure('Name',strrep(varName,'_',' '));
            hold;
            
            % collect the x values to plot each condition/sample at
            [xEggChamberVals, xEggChamberIDs] = ...
                obj.getSampleXValuesByEggChamber(eggChamberStagesToInclude);
            
            % generate color maps
            nc = numel(obj.condIndices);
            cmData =cbrewer('qual', 'Set1', max(nc,3));
            
            % collect the values of the metric for each condition/sample
            
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            for i=1:numel(obj.condIndices)
                xPlot = [];
                yPlot = [];
                s = obj.sampleIndices{i};
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude)
                            % collect values for the desired metric from all nuclei for the
                            % current sample/condition
                            x = obj.nucFullT.(varName)(...
                                obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                & obj.nucFullT.sample_Idx == s(j) ...
                                & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                    
                            % generate slightly offset x coordinates for each nucleus,
                            % centered around the sample X
                            nNuclei = size(x,1);
                            if nNuclei >1
                                % spacing between nuclei
                                nucSpacing = obj.spacingUnit ...
                                    * obj.cloudWidth / (nNuclei-1);
                    
                                % x coordinate for each nucleus of current condition/sample
                                curXPlot = xEggChamberVals{i}{j}(k) ...
                                    - floor(nNuclei/2)*nucSpacing ...
                                    + (0:(nNuclei-1))*nucSpacing;
                    
                                % y coordinate for each nucleus of current condition/sample
                                curYPlot = obj.nucFullT.(varName)(...
                                    obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucFullT.sample_Idx == s(j)...
                                    & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 1
                                curXPlot = xEggChamberVals{i}{j}(k);
                                curYPlot = obj.nucFullT.(varName)(...
                                    obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucFullT.sample_Idx == s(j)...
                                    & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 0
                                curXPlot = [];
                                curYPlot = [];
                            end
                    
                            % append coordinates of current condition/sample to global list
                            xPlot = [xPlot,curXPlot];
                            yPlot = [yPlot,curYPlot];

                            xErr = mean(curXPlot);
                            yErr = mean(curYPlot);
                            eErr = std(curYPlot);

                            xSampleValsVec = [xSampleValsVec,xEggChamberVals{i}{j}(k)];
                            xSampleIDsVec = [xSampleIDsVec,xEggChamberIDs{i}{j}{k}];

                            p = scatter(curXPlot,curYPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                            alpha(p,0.3);
                            errorbar(xErr,yErr,eErr,'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                        end
                    end
                end
                if ~isempty(yPlot)
                    yMax = max(yMax,max(yPlot(:)));
                    yMin = min(yMin,min(yPlot(:)));
                end
                
            end
            
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylim([min([0,1.2*yMin]),1.2*yMax]);
            ylabel(strrep(baseName,'_','\_'));
            grid on
        end

        %% scatter plot a metric by egg chamber (old)
        % older version, use scatterPlotAndSaveClusterMetricByEggChamber
        % instead
        function scatterPlotClustTableMetricByEggChamberOld(...
                obj,prefix,channel,baseName,suffix,eggChamberStagesToInclude,minClustVolume)
            % prefix: any allowable prefix which marks the compartment the metric is calculated on e.g. 'nuc', or 'clust'
            % channel: intensity channel , e.g. 1. 
                % (Value is ignored if the metric is a geometry feature 
                % rather than an intensity feature.
            % baseName: any allowable metric basename, e.g. 'Mean' or 'Volume'
            % suffix: any allowable suffix which marks processing steps,
            % e.g. 'raw' or 'eggChamberCorr'
            % if shorthand 'all' was used for cell stages to include, replace it by
            % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % check that the metric is present in the data table.
            varName = obj.buildVarName(prefix,channel,baseName,suffix,'geom');
            if ~ismember( varName, obj.clustT.Properties.VariableNames)
                disp(['Variable ',varName,' absent from table, cannot plot.']);
                return
            end

            % build figure
            figure('Name',strrep(varName,'_',' '));
            hold;
            
            % collect the x values to plot each condition/sample at
            [xEggChamberVals, xEggChamberIDs] = ...
                obj.getSampleXValuesByEggChamber(eggChamberStagesToInclude);
            
            % generate color maps
            nc = numel(obj.condIndices);
            cmData =cbrewer('qual', 'Set1', max(nc,3));
            
            % collect the values of the metric for each condition/sample
            
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            for i=1:numel(obj.condIndices)
                xPlot = [];
                yPlot = [];
                s = obj.sampleIndices{i};
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude)
                            % collect values for the desired metric from all nuclei for the
                            % current sample/condition
                            x = obj.clustT.(varName)(...
                                obj.clustT.cond_Idx ==obj.condIndices(i) ...
                                & obj.clustT.sample_Idx == s(j) ...
                                & obj.clustT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k) ...
                                & obj.clustT.clust_Volume >= minClustVolume);
                    
                            % generate slightly offset x coordinates for each nucleus,
                            % centered around the sample X
                            nClusters = size(x,1);
                            if nClusters >1
                                % spacing between nuclei
                                nucSpacing = obj.spacingUnit ...
                                    * obj.cloudWidth / (nClusters-1);
                    
                                % x coordinate for each nucleus of current condition/sample
                                curXPlot = xEggChamberVals{i}{j}(k) ...
                                    - floor(nClusters/2)*nucSpacing ...
                                    + (0:(nClusters-1))*nucSpacing;
                    
                                % y coordinate for each nucleus of current condition/sample
                                curYPlot = obj.clustT.(varName)(...
                                    obj.clustT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.clustT.sample_Idx == s(j)...
                                    & obj.clustT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k) ...
                                    & obj.clustT.clust_Volume >= minClustVolume)';
                    
                            elseif nClusters == 1
                                curXPlot = xEggChamberVals{i}{j}(k);
                                curYPlot = obj.clustT.(varName)(...
                                    obj.clustT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.clustT.sample_Idx == s(j)...
                                    & obj.clustT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k) ...
                                    & obj.clustT.clust_Volume >= minClustVolume)';
                    
                            elseif nClusters == 0
                                curXPlot = [];
                                curYPlot = [];
                            end
                    
                            % append coordinates of current condition/sample to global list
                            xPlot = [xPlot,curXPlot];
                            yPlot = [yPlot,curYPlot];

                            xErr = mean(curXPlot);
                            yErr = mean(curYPlot);
                            eErr = std(curYPlot);

                            xSampleValsVec = [xSampleValsVec,xEggChamberVals{i}{j}(k)];
                            xSampleIDsVec = [xSampleIDsVec,xEggChamberIDs{i}{j}{k}];

                            p = scatter(curXPlot,curYPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                            alpha(p,0.05);
                            errorbar(xErr,yErr,eErr,'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                        end
                    end
                end
                if ~isempty(yPlot)
                    yMax = max(yMax,max(yPlot(:)));
                    yMin = min(yMin,min(yPlot(:)));
                end
                
            end
            
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylim([min([0,1.2*yMin]),1.2*yMax]);
            ylabel(strrep(baseName,'_','\_'));
            grid on
        end
        
        %% scatter plot a nucleus metric by egg chamber
        function [nucTable,avgEcNucTable,avgCondNucTable,fh] = scatterPlotNucTableMetricByEggChamber(...
                obj,prefix,channel,baseName,suffix,metricName,...
                idxData,ecToExclude,conditionOrder,eggChamberStagesToInclude,removeBackgroundEggChambers,whichMetric,alphaVal)
            % prefix: any allowable prefix which marks the compartment the metric is calculated on e.g. 'nuc', or 'clust'
            % channel: intensity channel , e.g. 1. 
                % (Value is ignored if the metric is a geometry feature 
                % rather than an intensity feature.
            % baseName: any allowable metric basename, e.g. 'Mean' or 'Volume'
            % suffix: any allowable suffix which marks processing steps,
                % e.g. 'raw' or 'eggChamberCorr'
            % metricName: name used for the metric plotted in the display and table header
            % idxData: which clusters to plot. logical vector which size needs to match the the number of rows of
                % the nucTable
            % ecToExclude: egg chambers to exclude, formatted as {[i1,j1,k1],[i2,j2,k2],...}
                % where i is the index of the condition
                % j the index of the sample (FOV)
                % k the ID of the eggChamber
            % conditionOrder: array that sets the order of the clouds of points 
                % for the corresponding conditions from left to right, 
                % e.g. [1,3,2,4] plots condition index 1 first, then 3, then 2, then 4. 
            % eggChamberStagesToInclude: enter either 'all', or an array of
                % stages, e.g. [7,8] if you want to include only stages 7 and
                % 8.
            % minClustVolume: minimum Volume of Clusters (in um^3) to include in the plots
                % smaller clusters will be ignored
            % maxClustVolume: minimumVolume of Clusters (in um^3) to include in the plots
                % larger clusters will be ignored
            % removeBackgroundEggChambers: set to 1 to remove all egg
                % chambers with ID = 0; set to 0to include everything
            % whichMetric: set t0 'useMean' or 'useMedian' depending on
                % which summary metric you want to use.
            % alphaVal: value of the alpha (transparency) of individual
                % data points. 0 is fully transparent, 1 is fully solid.
            
            % check that the metric is present in the data table.
            varName = obj.buildVarName(prefix,channel,baseName,suffix,'geom');
            if ~ismember( varName, obj.nucFullT.Properties.VariableNames)
                disp(['Variable ',varName,' absent from table, cannot plot.']);
                return
            end

            [nucTable,avgEcNucTable,avgCondNucTable,fh] = ...
                obj.scatterPlotNucArbitraryMetricByEggChamber(...
                obj.nucFullT.(varName),metricName,idxData,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,...
                removeBackgroundEggChambers,whichMetric,alphaVal);
        end

        %% scatter plot any nucleus metric by egg chamber
        function [nucTable,avgEcNucTable,avgCondNucTable,fh] = scatterPlotNucArbitraryMetricByEggChamber(...
                obj,yData,metricName,idxData,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,removeBackgroundEggChambers,whichMetric,alphaVal)
            
            % yData: the data to plot. Vector of values which size needs to match the number of rows of
            % the nucTable
            % idxData: which clusters to plot. logical vector which size needs to match the the number of rows of
            % the nucTable
            % metricName: name used for the metric plotted in the display and table header
            % ecToExclude: egg chambers to exclude, formatted as {[i1,j1,k1],[i2,j2,k2],...}
                % where i is the index of the condition
                % j the index of the sample (FOV)
                % k the ID of the eggChamber
            % conditionOrder: array that sets the order of the clouds of points 
                % for the corresponding conditions from left to right, 
                % e.g. [1,3,2,4] plots condition index 1 first, then 3, then 2, then 4. 
            % eggChamberStagesToInclude: enter either 'all', or an array of
                % stages, e.g. [7,8] if you want to include only stages 7 and
                % 8.
            % minClustVolume: minimum Volume of Clusters (in um^3) to include in the plots
                % smaller clusters will be ignored
            % maxClustVolume: minimumVolume of Clusters (in um^3) to include in the plots
                % larger clusters will be ignored
            % removeBackgroundEggChambers: set to 1 to remove all egg
                % chambers with ID = 0; set to 0to include everything
            % whichMetric: set t0 'useMean' or 'useMedian' depending on
                % which summary metric you want to use.
            % alphaVal: value of the alpha (transparency) of individual
                % data points. 0 is fully transparent, 1 is fully solid.


            % if shorthand 'all' was used for cell stages to include, replace it by
                % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end

            % reformat the indices to exclude from cell {[i1,j1,k1], [i2,j2,k2],...}
            % to matrix    [i1,j1,k1; 
            %               i2,j2,k2; ...]
            
            if isempty(ecToExclude)
                idxToExclude = [NaN,NaN,NaN];
            else
                idxToExclude = zeros(numel(ecToExclude),3);
                for i=1:numel(ecToExclude)
                    if numel(ecToExclude{i})~=3
                        disp(["ecToExclude entry #",num2str(i)," contains the wrong number of elements (",...
                            nu2mstr(numel(ecToExclude{i})),"); it should have three elements: ",...
                            "first the condition to exclude, second the FOV to exclude,",...
                            " third the egg chamber to exclude"]);
                        return
                    else
                        idxToExclude(i,:) = ecToExclude{i};
                    end 
                end
            end

            % build figure
            if ~isempty(metricName)
                figName = metricName;
            else
                figName = strrep(varName,'_',' ');
            end
            fh = figure('Name',figName);
            hold;
             
            % generate color maps
            nc = numel(obj.condIndices);
            cmData =cbrewer('qual', 'Dark2', max(nc,3));
            if nc ==2
                cmData = cmData([1,3],:);
            end

            % collect the values of the metric for each condition/sample
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            xCtr = 0;
            curCond = 0;

            % initialize variables that will be output as tables
            % full table (all nuclei)
            nTableCondition_Name = {};
            nTableSample = [];
            nTableFileName = {};
            nTableEggChamberID = [];
            nTableEggChamberStage = [];
            nTableNucleusLabel = [];
            nTableXData = [];
            nTableYData = [];

             % average over each egg chamber
            avgEcTableCondition_Name = {};
            avgEcTableSample = [];
            avgEcTableFileName = {};
            avgEcTableEggChamberID = [];
            avgEcTableEggChamberStage = [];
            avgEcTableXData = [];
            avgEcTableYData = [];
            avgEcTableSTDData = [];
            avgEcTableSEMData = [];
            avgEcTableNData = [];
            
            % average over each condition
            avgCondTableCondition = {};
            avgCondTableXData = [];
            avgCondTableYData = [];
            avgCondTableSTDData = [];
            avgCondTableSEMData = [];
            avgCondTableNData = [];
            for idxCondOrder=1:numel(conditionOrder)
                i = find(ismember(obj.condIndices,conditionOrder(idxCondOrder)));
                if ~isempty(i)
                    xPlot = [];
                    yPlot = [];
                    s = obj.sampleIndices{i};
                    for j=1:obj.nSamples(i)
                        for k=1:obj.eggChamberNumber{i}(j)
                            % figure out if current eggchamber needs to be
                            % filtered out
                            if removeBackgroundEggChambers == 1
                                processEggChamber = ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude) ...
                                    && ~ismember([obj.condIndices(i),obj.sampleIndices{i}(j),obj.eggChamberIDs{i}{j}(k)],...
                                    idxToExclude,'rows')...
                                    && obj.eggChamberIDs{i}{j}(k)~=0;  
                            else
                                processEggChamber = ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude) ...
                                    && ~ismember([obj.condIndices(i),obj.sampleIndices{i}(j),obj.eggChamberIDs{i}{j}(k)],...
                                    idxToExclude,'rows');
                            end
                            if processEggChamber
    
                                % collect values for the desired metric from all nuclei for the
                                % current sample/condition
                                x = yData(idxData ...
                                    & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucFullT.sample_Idx == s(j) ...
                                    & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                        
                                % generate slightly offset x coordinates for each nucleus,
                                % centered around the sample X
                                if (curCond ~= 0) && (curCond ~= obj.condIndices(i)) % add extra offset for each condition change
                                    % plot average/std of previous
                                    % condition
                                    if sum(~isnan(yCondAvg)) > 0
                                        if sum(~isnan(yCondAvg)) > 1
                                            nDenom = sqrt(sum(~isnan(yCondAvg) & ~isinf(yCondAvg))-1);
                                        else
                                            nDenom = 1;
                                        end
                                        if strcmp(whichMetric,'useMean')
                                            errorbar(mean(xCondAvg,'omitnan'),mean(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                                                std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom,...
                                                'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                                            avgCondTableYData = [avgCondTableYData;...
                                                mean(yCondAvg(~isinf(yCondAvg)),'omitnan')];

                                        elseif strcmp(whichMetric,'useMedian')
                                            errorbar(mean(xCondAvg,'omitnan'),median(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                                                std(yCondAvg,'omitnan')/nDenom,...
                                                'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                                            avgCondTableYData = [avgCondTableYData;...
                                                median(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                                        else
                                            disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                                        end

                                        avgCondTableCondition = [avgCondTableCondition;...
                                            obj.conditionNames(find(obj.condIndices == curCond,1))];
                                        avgCondTableXData = [avgCondTableXData;...
                                                mean(xCondAvg,'omitnan')];
                                        avgCondTableSTDData = [avgCondTableSTDData;...
                                            std(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                                        avgCondTableSEMData = [avgCondTableSEMData;...
                                            std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom];
                                        avgCondTableNData = [avgCondTableNData;...
                                            sum(~isnan(yCondAvg) & ~isinf(yCondAvg))];
                                        
                                    else
                                        avgCondTableCondition = [avgCondTableCondition;...
                                            obj.conditionNames(find(obj.condIndices == curCond,1 ))];
                                        avgCondTableXData = [avgCondTableXData;NaN];
                                        avgCondTableYData = [avgCondTableYData;NaN];
                                        avgCondTableSTDData = [avgCondTableSTDData;NaN];
                                        avgCondTableSEMData = [avgCondTableSEMData;NaN];
                                        avgCondTableNData = [avgCondTableNData;0];
                                        
                                    end

                                    curCond = obj.condIndices(i);
                                    xCtr = xCtr + obj.condSeparator*obj.spacingUnit;
                                    xCondAvg = xCtr;
                                    yCondAvg = [];
                                elseif curCond == 0
                                    curCond = obj.condIndices(i);
                                    xCondAvg = xCtr;
                                    yCondAvg = [];
                                else
                                    xCtr = xCtr+obj.ecSeparator*obj.spacingUnit;
                                    xCondAvg = [xCondAvg;xCtr];
                                end
                                
                                nNuc = size(x,1);
                                if nNuc >1
                                    nucSpacing = obj.cloudWidth*obj.spacingUnit / (nNuc-1);
                        
                                    % x coordinate for each nucleus of current condition/sample
                                    curXPlot = xCtr ...
                                        - floor(nNuc/2)*nucSpacing ...
                                        + (0:(nNuc-1))*nucSpacing;
                        
                                    % y coordinate for each nucleus of current condition/sample
                                    curYPlot = yData(idxData ...
                                        & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                        & obj.nucFullT.sample_Idx == s(j)...
                                        & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                        
                                elseif nNuc == 1
                                    curXPlot = xCtr;
                                    curYPlot = yData(idxData ...
                                        & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                        & obj.nucFullT.sample_Idx == s(j)...
                                        & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                        
                                elseif nNuc == 0
                                    curXPlot = [];
                                    curYPlot = [];
                                end
                                
                                % fill in cluster Table values
                                nTableCondition_Name = [nTableCondition_Name; ...
                                    repmat(obj.conditionNames(i),size(curXPlot'))];

                                nTableSample = [nTableSample; ...
                                    repmat(s(j),size(curXPlot'))];

                                curFname = obj.nucFullT.sample_InputFileName(idxData ...
                                        & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                        & obj.nucFullT.sample_Idx == s(j)...
                                        & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                                nTableFileName = [nTableFileName; curFname];

                                nTableEggChamberID = [nTableEggChamberID; ...
                                    repmat(obj.eggChamberIDs{i}{j}(k),size(curXPlot'))];

                                curStage = obj.nucFullT.eggChamber_Stage(idxData ...
                                        & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                        & obj.nucFullT.sample_Idx == s(j)...
                                        & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                                nTableEggChamberStage = [nTableEggChamberStage;curStage];
                                
                                curNuc = obj.nucFullT.nuc_Label(idxData ...
                                        & obj.nucFullT.cond_Idx ==obj.condIndices(i) ...
                                        & obj.nucFullT.sample_Idx == s(j)...
                                        & obj.nucFullT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                                nTableNucleusLabel = [nTableNucleusLabel; curNuc];
                                
                                nTableXData = [nTableXData;curXPlot'];
                                nTableYData = [nTableYData;curYPlot'];
                        
                                % append coordinates of current condition/sample to global list
                                xPlot = [xPlot,curXPlot];
                                yPlot = [yPlot,curYPlot];
                                
                                xMeanEC = mean(curXPlot);
                                if strcmp(whichMetric,'useMean')
                                    yMeanEC = mean(curYPlot(~isinf(curYPlot)),'omitnan');
                                elseif strcmp(whichMetric,'useMedian')
                                    yMeanEC = median(curYPlot(~isinf(curYPlot)),'omitnan');
                                else
                                    disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                                end
                                
                                yCondAvg = [yCondAvg;yMeanEC];
                                curID = [obj.conditionNames{i},...
                                    ' sample\_',num2str(obj.sampleIndices{i}(j)),...
                                    ' eggChamber\_',num2str(obj.eggChamberIDs{i}{j}(k))];
                                xSampleValsVec = [xSampleValsVec,xCtr];
                                xSampleIDsVec = [xSampleIDsVec,curID];

                                avgEcTableCondition_Name = [avgEcTableCondition_Name; obj.conditionNames{i}];
                                avgEcTableSample = [avgEcTableSample; obj.sampleIndices{i}(j)];
                                avgEcTableFileName = [avgEcTableFileName; curFname{1}];
                                avgEcTableEggChamberID = [avgEcTableEggChamberID; obj.eggChamberIDs{i}{j}(k)];
                                avgEcTableEggChamberStage = [avgEcTableEggChamberStage; obj.eggChamberStages{i}{j}(k)];
                                avgEcTableXData = [avgEcTableXData;xMeanEC];
                                avgEcTableYData = [avgEcTableYData; yMeanEC];
                                avgEcTableSTDData = [avgEcTableSTDData; std(curYPlot)];
                                if numel(curXPlot) == 0
                                    avgEcTableSEMData = [avgEcTableSEMData; NaN];
                                elseif numel(curXPlot) == 1
                                    avgEcTableSEMData = [avgEcTableSEMData; std(curYPlot)];
                                else
                                    avgEcTableSEMData = [avgEcTableSEMData; std(curYPlot)/sqrt(numel(curXPlot)-1)];
                                end
                                avgEcTableNData = [avgEcTableNData;numel(curXPlot)];
    
                                p = scatter(curXPlot,curYPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                                alpha(p,alphaVal);
                                p2 = scatter(xMeanEC,yMeanEC,'d','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                                alpha(p2,0.9);
                                %errorbar(xErr,yErr,eErr,'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                            end
                        end
                    end
                    if ~isempty(yPlot)
                        yMax = max(yMax,max(yPlot(:)));
                        yMin = min(yMin,min(yPlot(:)));
                    end
                end
            end

            % plot average of last condition
            if sum(~isnan(yCondAvg) & ~isinf(yCondAvg)) > 0
                if sum(~isnan(yCondAvg) & ~isinf(yCondAvg)) > 1
                    nDenom = sqrt(sum(~isnan(yCondAvg) & ~isinf(yCondAvg))-1);
                else
                    nDenom = 1;
                end
                if strcmp(whichMetric,'useMean')
                    errorbar(mean(xCondAvg),mean(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                        std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom,...
                        'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                    avgCondTableYData = [avgCondTableYData;...
                        mean(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                elseif strcmp(whichMetric,'useMedian')
                    errorbar(mean(xCondAvg),median(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                        std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom,...
                        'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                    avgCondTableYData = [avgCondTableYData;...
                        median(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                else
                    disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                end
                avgCondTableCondition = [avgCondTableCondition;...
                    obj.conditionNames(find(obj.condIndices == curCond,1))];
                avgCondTableXData = [avgCondTableXData;...
                    mean(xCondAvg)];
                avgCondTableSTDData = [avgCondTableSTDData;...
                    std(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                avgCondTableSEMData = [avgCondTableSEMData;...
                    std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom];
                avgCondTableNData = [avgCondTableNData;...
                    sum(~isnan(yCondAvg(~isinf(yCondAvg))))];
            else
                avgCondTableCondition = [avgCondTableCondition;...
                    obj.conditionNames(find(obj.condIndices == curCond,1))];
                avgCondTableXData = [avgCondTableXData;NaN];
                avgCondTableYData = [avgCondTableYData;NaN];
                avgCondTableSTDData = [avgCondTableSTDData;NaN];
                avgCondTableSEMData = [avgCondTableSEMData;NaN];
                avgCondTableNData = [avgCondTableNData;0];
            end
                                    
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);

            if isinf(yMin) && yMin >= 0
                if isinf(yMax) && yMax <= 0
                    yMin = 0;
                    yMax = 1;
                end
            else
                yMin = min([0,1.2*yMin]);
                yMax = 1.2*yMax;
                if yMin == yMax
                    if yMin == 0
                        yMin = yMin - 1;
                        yMax = yMax + 1;
                    else
                        yMin = yMin - 0.1*yMin;
                        yMax = yMax + 0.1*yMin;
                    end
                end
            end  
            ylim([yMin,yMax]);
            
            ylabel(strrep(metricName,'_','\_'));
            grid on

            % build full table
            nucTable = table(nTableCondition_Name,nTableSample,nTableFileName,nTableEggChamberID,...
                nTableEggChamberStage,nTableNucleusLabel, ...
                nTableXData,nTableYData,'VariableNames',{...
                'Condition','sample','sample_InputFileName','eggChamber','eggChamber_Stage',...
                'nuc_Label','x',metricName});  

            % build table of averages by eggChamber
            avgEcNucTable = table(avgEcTableCondition_Name, avgEcTableSample, avgEcTableFileName,...
                avgEcTableEggChamberID, avgEcTableEggChamberStage, avgEcTableXData, avgEcTableYData, ...
                avgEcTableSTDData, avgEcTableSEMData,avgEcTableNData,'VariableNames',...
                {'Condition','sample','sample_InputFileName','eggChamber','eggChamber_Stage',...
                'x',[metricName,'_EggChamberMean'],[metricName,'_EggChamberSTD'],[metricName,'_EggChamberSEM'],...
                'N_nuclei'});

            % build average table
            avgCondNucTable = table(avgCondTableCondition, avgCondTableXData, avgCondTableYData, ...
                avgCondTableSTDData, avgCondTableSEMData,avgCondTableNData,'VariableNames',...
                {'Condition','x',[metricName,'_Mean'],[metricName,'_STD'],[metricName,'_SEM'],...
                'N_eggChambers'});
        end

        %% scatter plot a metric from the cluster table by egg chamber
        function [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fh] = scatterPlotClustTableMetricByEggChamber(...
                obj,prefix,channel,baseName,suffix,metricName,...
                idxData,ecToExclude,conditionOrder,eggChamberStagesToInclude,...
                minClustVolume,maxClustVolume,removeBackgroundEggChambers,whichMetric,alphaVal)
            % prefix: any allowable prefix which marks the compartment the metric is calculated on e.g. 'nuc', or 'clust'
            % channel: intensity channel , e.g. 1. 
                % (Value is ignored if the metric is a geometry feature 
                % rather than an intensity feature.
            % baseName: any allowable metric basename, e.g. 'Mean' or 'Volume'
            % suffix: any allowable suffix which marks processing steps,
                % e.g. 'raw' or 'eggChamberCorr'
            % metricName: name used for the metric plotted in the display and table header
            % ecToExclude: egg chambers to exclude, formatted as {[i1,j1,k1],[i2,j2,k2],...}
                % where i is the index of the condition
                % j the index of the sample (FOV)
                % k the ID of the eggChamber
            % conditionOrder: array that sets the order of the clouds of points 
                % for the corresponding conditions from left to right, 
                % e.g. [1,3,2,4] plots condition index 1 first, then 3, then 2, then 4. 
            % eggChamberStagesToInclude: enter either 'all', or an array of
                % stages, e.g. [7,8] if you want to include only stages 7 and
                % 8.
            % minClustVolume: minimum Volume of Clusters (in um^3) to include in the plots
                % smaller clusters will be ignored
            % maxClustVolume: minimumVolume of Clusters (in um^3) to include in the plots
                % larger clusters will be ignored
            % removeBackgroundEggChambers: set to 1 to remove all egg
                % chambers with ID = 0; set to 0to include everything
            % whichMetric: set t0 'useMean' or 'useMedian' depending on
                % which summary metric you want to use.
            % alphaVal: value of the alpha (transparency) of individual
                % data points. 0 is fully transparent, 1 is fully solid.
            
            % check that the metric is present in the data table.
            varName = obj.buildVarName(prefix,channel,baseName,suffix,'geom');
            if ~ismember( varName, obj.clustT.Properties.VariableNames)
                disp(['Variable ',varName,' absent from table, cannot plot.']);
                return
            end

            [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fh] = ...
                obj.scatterPlotClustArbitraryMetricByEggChamber(obj.clustT.(varName),metricName,idxData,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,...
                minClustVolume,maxClustVolume,removeBackgroundEggChambers,whichMetric,alphaVal);
        end

        %% scatter plot a metric from the cluster table by egg chamber, stratified by index
        function [clustTable1,avgNucClustTable1,avgEcClustTable1,avgCondClustTable1,...
                clustTable2,avgNucClustTable2,avgEcClustTable2,avgCondClustTable2,fh] = ...
                scatterPlotClustDualArbitraryMetricByEggChamber(...
                obj,yData,metricName,idxData1,idxData2,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,...
                minClustVolume,maxClustVolume,removeBackgroundEggChambers,whichMetric,alphaVal,userDefinedColorMap)
            % prefix: any allowable prefix which marks the compartment the metric is calculated on e.g. 'nuc', or 'clust'
            % channel: intensity channel , e.g. 1. 
                % (Value is ignored if the metric is a geometry feature 
                % rather than an intensity feature.
            % baseName: any allowable metric basename, e.g. 'Mean' or 'Volume'
            % suffix: any allowable suffix which marks processing steps,
                % e.g. 'raw' or 'eggChamberCorr'
            % metricName: name used for the metric plotted in the display and table header
            % ecToExclude: egg chambers to exclude, formatted as {[i1,j1,k1],[i2,j2,k2],...}
                % where i is the index of the condition
                % j the index of the sample (FOV)
                % k the ID of the eggChamber
            % conditionOrder: array that sets the order of the clouds of points 
                % for the corresponding conditions from left to right, 
                % e.g. [1,3,2,4] plots condition index 1 first, then 3, then 2, then 4. 
            % eggChamberStagesToInclude: enter either 'all', or an array of
                % stages, e.g. [7,8] if you want to include only stages 7 and
                % 8.
            % minClustVolume: minimum Volume of Clusters (in um^3) to include in the plots
                % smaller clusters will be ignored
            % maxClustVolume: minimumVolume of Clusters (in um^3) to include in the plots
                % larger clusters will be ignored
            % removeBackgroundEggChambers: set to 1 to remove all egg
                % chambers with ID = 0; set to 0to include everything
            % whichMetric: set t0 'useMean' or 'useMedian' depending on
                % which summary metric you want to use.
            % alphaVal: value of the alpha (transparency) of individual
                % data points. 0 is fully transparent, 1 is fully solid.
            % userdefinedColorMap: you can enter a color value - this will supersede the autocolors FOR ALL. 
                % Enter a 2x1 cell array, where each cell is an [nConditions x 3 matrix].         

            obj.setPlottingDistances(2,3,1,1,0.05);

            xOffset1 = -obj.ecSeparator/4;   
            [clustTable1,avgNucClustTable1,avgEcClustTable1,avgCondClustTable1,fh] = ...
                scatterPlotClustArbitraryMetricByEggChamber(obj,yData,metricName,idxData1,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,minClustVolume,maxClustVolume,...
                removeBackgroundEggChambers,whichMetric,alphaVal,[],xOffset1,userDefinedColorMap{1});
            
            xOffset2 = obj.ecSeparator/4;
            [clustTable2,avgNucClustTable2,avgEcClustTable2,avgCondClustTable2,fh] = ...
                scatterPlotClustArbitraryMetricByEggChamber(obj,yData,metricName,idxData2,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,minClustVolume,maxClustVolume,...
                removeBackgroundEggChambers,whichMetric,alphaVal,fh,xOffset2,userDefinedColorMap{2});
            
            obj.resetPlottingDistances();
        end
        
        %% scatter plot any cluster metric by egg chamber
        function [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fh] = scatterPlotClustArbitraryMetricByEggChamber(...
                obj,yData,metricName,idxData,...
                ecToExclude,conditionOrder,eggChamberStagesToInclude,minClustVolume,maxClustVolume,...
                removeBackgroundEggChambers,whichMetric,alphaVal,varargin)
            % yData: the data to plot. Vector of values which size needs to match the number of rows of
            % the clustTable
            % idxData: which clusters to plot. logical vector which size needs to match the the number of rows of
            % the clustTable
            % metricName: name used for the metric plotted in the display and table header
            % ecToExclude: egg chambers to exclude, formatted as {[i1,j1,k1],[i2,j2,k2],...}
                % where i is the index of the condition
                % j the index of the sample (FOV)
                % k the ID of the eggChamber
            % conditionOrder: array that sets the order of the clouds of points 
                % for the corresponding conditions from left to right, 
                % e.g. [1,3,2,4] plots condition index 1 first, then 3, then 2, then 4. 
            % eggChamberStagesToInclude: enter either 'all', or an array of
                % stages, e.g. [7,8] if you want to include only stages 7 and
                % 8.
            % minClustVolume: minimum Volume of Clusters (in um^3) to include in the plots
                % smaller clusters will be ignored
            % maxClustVolume: minimumVolume of Clusters (in um^3) to include in the plots
                % larger clusters will be ignored
            % removeBackgroundEggChambers: set to 1 to remove all egg
                % chambers with ID = 0; set to 0to include everything
            % whichMetric: set t0 'useMean' or 'useMedian' depending on
                % which summary metric you want to use.
            % alphaVal: value of the alpha (transparency) of individual
                % data points. 0 is fully transparent, 1 is fully solid.
            % optional arguments varargin: 
            % #1 you can enter the handle of a
                % figure to overlay the new data onto it
            % #2 you can enter an x offset value - this will translate all
                % the data points by the amount specified by offset (used when plotting two subconditions per EC).
            % #3 you can enter a color value - this will supersede the autocolors FOR ALL. 
                % Enter a nConditions x 3 matrix.     

            %collect figure handle if present   
            userDefinedColorMap = [];
            if ~isempty(varargin)
                fh =varargin{1};
                if numel(varargin)>1
                    xOffset = varargin{2};
                    if numel(varargin)>2
                        userDefinedColorMap = varargin{3};
                    end
                end
            else
                fh = [];
                xOffset = 0;
            end

            % if shorthand 'all' was used for cell stages to include, replace it by
                % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % check that the metric and index have the right size.
            if numel(yData) ~= size(obj.clustT,1)
                disp(['yData has ',num2str(numel(yData)),...
                    ' elements, which does not match the size of the cluster table (',...
                    num2str(size(obj.clustT,1)),')']);
                return
            end
            if numel(idxData) ~= size(obj.clustT,1)
                disp(['idxData has ',num2str(numel(idxData)),...
                    ' elements, which does not match the size of the cluster table (',...
                    num2str(size(obj.clustT,1)),')']);
                return
            end

            % reformat the indices to exclude from cell {[i1,j1,k1], [i2,j2,k2],...}
            % to matrix    [i1,j1,k1; 
            %               i2,j2,k2; ...]
            
            if isempty(ecToExclude)
                idxToExclude = [NaN,NaN,NaN];
            else
                idxToExclude = zeros(numel(ecToExclude),3);
                for i=1:numel(ecToExclude)
                    if numel(ecToExclude{i})~=3
                        disp(["ecToExclude entry #",num2str(i)," contains the wrong number of elements (",...
                            nu2mstr(numel(ecToExclude{i})),"); it should have three elements: ",...
                            "first the condition to exclude, second the FOV to exclude,",...
                            " third the egg chamber to exclude"]);
                        return
                    else
                        idxToExclude(i,:) = ecToExclude{i};
                    end 
                end
            end

            % build figure
            if ~isempty(metricName)
                figName = [metricName,' by Clusters'];
            else
                figName = [strrep(varName,'_',' '),' by Clusters'];
            end

            if isempty(fh)
                fh = figure('Name',figName);
                hold;
                 % generate color maps
                nc = numel(obj.condIndices);
                if isempty(userDefinedColorMap)
                    cmData =cbrewer('qual', 'Dark2', max(nc,3));
                    if nc ==2
                        cmData = cmData([1,3],:);
                    end
                else
                    cmData = userDefinedColorMap;
                end 
            else
                 % generate color maps
                nc = numel(obj.condIndices);
                if isempty(userDefinedColorMap)
                    cmData =cbrewer('qual', 'Dark2', max(nc,3));
                    if nc ==2
                        cmData = cmData([1,3],:);
                    end
                else
                    cmData = userDefinedColorMap;
                end
                
            end
           
            % collect the values of the metric for each condition/sample
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            xCtr = xOffset;
            curCond = 0;

            % initialize variables that will be output as tables
            % full table (all data points)
            cTableCondition_Name = {};
            cTableSample = [];
            cTableFileName = {};
            cTableEggChamberID = [];
            cTableEggChamberStage = [];
            cTableNucleusLabel = [];
            cTableCLustLabel = [];
            cTableXData = [];
            cTableYData = [];
            
            % average over each egg chamber
            avgEcTableCondition_Name = {};
            avgEcTableSample = [];
            avgEcTableFileName = {};
            avgEcTableEggChamberID = [];
            avgEcTableEggChamberStage = [];
            avgEcTableXData = [];
            avgEcTableYData = [];
            avgEcTableSTDData = [];
            avgEcTableSEMData = [];
            avgEcTableNData = [];
            
            % average over each condition
            avgCondTableCondition_Name = {};
            avgCondTableXData = [];
            avgCondTableYData = [];
            avgCondTableSTDData = [];
            avgCondTableSEMData = [];
            avgCondTableNData = [];
            for idxCondOrder=1:numel(conditionOrder)
                i = find(ismember(obj.condIndices,conditionOrder(idxCondOrder)));
                if ~isempty(i)
                    xPlot = [];
                    yPlot = [];
                    s = obj.sampleIndices{i};
                    for j=1:obj.nSamples(i)
                        for k=1:obj.eggChamberNumber{i}(j)
                            % figure out if current eggchamber needs to be
                            % filtered out
                            if removeBackgroundEggChambers == 1
                                processEggChamber = ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude) ...
                                    && ~ismember([obj.condIndices(i),obj.sampleIndices{i}(j),obj.eggChamberIDs{i}{j}(k)],...
                                    idxToExclude,'rows')...
                                    && obj.eggChamberIDs{i}{j}(k)~=0;  
                            else
                                processEggChamber = ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude) ...
                                    && ~ismember([obj.condIndices(i),obj.sampleIndices{i}(j),obj.eggChamberIDs{i}{j}(k)],...
                                    idxToExclude,'rows');
                            end
                            if processEggChamber
                                % collect values for the desired metric from all nuclei for the
                                % current sample/condition
                                curIdx = idxData ...
                                    & obj.clustT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.clustT.sample_Idx == s(j) ...
                                    & obj.clustT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k) ...
                                    & obj.clustT.clust_Volume >= minClustVolume...
                                    & obj.clustT.clust_Volume <= maxClustVolume;

                                x = yData(curIdx);
                        
                                % generate slightly offset x coordinates for each nucleus,
                                % centered around the sample X
                                if (curCond ~= 0) && (curCond ~= obj.condIndices(i)) % add extra offset for each condition change
                                    % plot average/std of previous
                                    % condition
                                    if sum(~isnan(yCondAvg)) > 0
                                        if sum(~isnan(yCondAvg)) > 1
                                            nDenom = sqrt(sum(~isnan(yCondAvg))-1);
                                        else
                                            nDenom = 1;
                                        end
                                        if strcmp(whichMetric,'useMean')
                                            errorbar(mean(xCondAvg,'omitnan'),mean(yCondAvg,'omitnan'),std(yCondAvg,'omitnan')/nDenom,...
                                                'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                                            avgCondTableYData = [avgCondTableYData;...
                                                mean(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                                        elseif strcmp(whichMetric,'useMedian')
                                            errorbar(mean(xCondAvg,'omitnan'),median(yCondAvg,'omitnan'),std(yCondAvg,'omitnan')/nDenom,...
                                                'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                                            avgCondTableYData = [avgCondTableYData;...
                                                median(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                                        else
                                            disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                                        end
                                        avgCondTableCondition_Name = [avgCondTableCondition_Name;...
                                            obj.conditionNames(find(obj.condIndices == curCond,1))];
                                        avgCondTableXData = [avgCondTableXData;...
                                            mean(xCondAvg,'omitnan')];
                                        avgCondTableSTDData = [avgCondTableSTDData;...
                                            std(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                                        avgCondTableSEMData = [avgCondTableSEMData;...
                                            std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom];
                                        avgCondTableNData = [avgCondTableNData;...
                                            sum(~isnan(yCondAvg) & ~isinf(yCondAvg))];
                                        
                                    else
                                        avgCondTableCondition_Name = [avgCondTableCondition_Name;...
                                            obj.conditionNames(find(obj.condIndices == curCond,1 ))];
                                        avgCondTableXData = [avgCondTableXData;NaN];
                                        avgCondTableYData = [avgCondTableYData;NaN];
                                        avgCondTableSTDData = [avgCondTableSTDData;NaN];
                                        avgCondTableSEMData = [avgCondTableSEMData;NaN];
                                        avgCondTableNData = [avgCondTableNData;0];             
                                    end

                                    curCond = obj.condIndices(i);
                                    xCtr = xCtr+obj.condSeparator*obj.spacingUnit;
                                    xCondAvg = xCtr;
                                    yCondAvg = [];
                                elseif curCond == 0
                                    curCond = obj.condIndices(i);
                                    xCondAvg = xCtr;
                                    yCondAvg = [];
                                else
                                    xCtr = xCtr+obj.ecSeparator;
                                    xCondAvg = [xCondAvg;xCtr];
                                end
                                
                                nClusters = size(x,1);
                                if nClusters >1
                                    % spacing between nuclei
                                    nucSpacing = obj.cloudWidth*obj.spacingUnit / (nClusters-1);
                                    
                                    % x coordinate for each nucleus of current condition/sample
                                    curXPlot = xCtr ...
                                        - floor(nClusters/2)*nucSpacing ...
                                        + (0:(nClusters-1))*nucSpacing;
                        
                                    % y coordinate for each nucleus of current condition/sample
                                    curYPlot = yData(curIdx)';
                        
                                elseif nClusters == 1
                                    curXPlot = xCtr;
                                    curYPlot = yData(curIdx)';
                        
                                elseif nClusters == 0
                                    curXPlot = [];
                                    curYPlot = [];
                                end
                                
                                % fill in cluster Table values
                                cTableCondition_Name = [cTableCondition_Name; ...
                                    repmat(obj.conditionNames(i),size(curXPlot'))];

                                cTableSample = [cTableSample; ...
                                    repmat(s(j),size(curXPlot'))];
                                
                                curFname = obj.clustT.sample_InputFileName(curIdx);
                                cTableFileName = [cTableFileName; curFname];

                                cTableEggChamberID = [cTableEggChamberID; ...
                                    repmat(obj.eggChamberIDs{i}{j}(k),size(curXPlot'))];

                                curStage = obj.clustT.eggChamber_Stage(curIdx);
                                cTableEggChamberStage = [cTableEggChamberStage;curStage];
                                
                                curNuc = obj.clustT.nuc_Label(curIdx);
                                cTableNucleusLabel = [cTableNucleusLabel; curNuc];
                                
                                curClust = obj.clustT.clust_Label(curIdx);
                                cTableCLustLabel = [cTableCLustLabel; curClust];

                                cTableXData = [cTableXData;curXPlot'];
                                cTableYData = [cTableYData;curYPlot'];
                        
                                % append coordinates of current condition/sample to global list
                                xPlot = [xPlot,curXPlot];
                                yPlot = [yPlot,curYPlot];
                                
                                xMeanEC = mean(curXPlot);
                                if strcmp(whichMetric,'useMean')
                                    yMeanEC = mean(curYPlot(~isinf(curYPlot)),'omitnan');
                                elseif strcmp(whichMetric,'useMedian')
                                    yMeanEC = median(curYPlot(~isinf(curYPlot)),'omitnan');
                                else
                                    disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                                end

                                yCondAvg = [yCondAvg;yMeanEC];
                                curID = [obj.conditionNames{i},...
                                    ' sample\_',num2str(obj.sampleIndices{i}(j)),...
                                    ' eggChamber\_',num2str(obj.eggChamberIDs{i}{j}(k))];
                                xSampleValsVec = [xSampleValsVec,xCtr];
                                xSampleIDsVec = [xSampleIDsVec,curID];
                                
                                avgEcTableCondition_Name = [avgEcTableCondition_Name; obj.conditionNames{i}];
                                avgEcTableSample = [avgEcTableSample; obj.sampleIndices{i}(j)];
                                if numel(curFname)>=1
                                    avgEcTableFileName = [avgEcTableFileName; curFname{1}];
                                else
                                    x = obj.clustT.sample_InputFileName(obj.clustT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.clustT.sample_Idx == s(j));
                                    if ~isempty(x)
                                        avgEcTableFileName = [avgEcTableFileName; x{1}];
                                    else
                                        avgEcTableFileName = [avgEcTableFileName; 'unknown'];
                                    end
                                    
                                end
                                
                                avgEcTableEggChamberID = [avgEcTableEggChamberID; obj.eggChamberIDs{i}{j}(k)];
                                avgEcTableEggChamberStage = [avgEcTableEggChamberStage; obj.eggChamberStages{i}{j}(k)];
                                avgEcTableXData = [avgEcTableXData;xMeanEC];
                                avgEcTableYData = [avgEcTableYData; yMeanEC];
                                avgEcTableSTDData = [avgEcTableSTDData; std(curYPlot)];
                                if numel(curXPlot) == 0
                                    avgEcTableSEMData = [avgEcTableSEMData; NaN];
                                elseif numel(curXPlot) == 1
                                    avgEcTableSEMData = [avgEcTableSEMData; std(curYPlot)];
                                else
                                    avgEcTableSEMData = [avgEcTableSEMData; std(curYPlot)/sqrt(numel(curXPlot)-1)];
                                end
                                avgEcTableNData = [avgEcTableNData;numel(curXPlot)];

                                p = scatter(curXPlot,curYPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                                alpha(p,alphaVal);
                                p2 = scatter(xMeanEC,yMeanEC,'d','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                                alpha(p2,0.9);
                                %errorbar(xErr,yErr,eErr,'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                            end
                        end
                    end
                    if ~isempty(yPlot)
                        yMax = max(yMax,max(yPlot(:)));
                        yMin = min(yMin,min(yPlot(:)));
                    end
                end
            end

            % plot average of last condition
            if sum(~isnan(yCondAvg) & ~isinf(yCondAvg)) > 0
                if sum(~isnan(yCondAvg) & ~isinf(yCondAvg)) > 1
                    nDenom = sqrt(sum(~isnan(yCondAvg) & ~isinf(yCondAvg))-1);
                else
                    nDenom = 1;
                end
                if strcmp(whichMetric,'useMean')
                    errorbar(mean(xCondAvg),mean(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                        std(yCondAvg,'omitnan')/nDenom,...
                        'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                    avgCondTableYData = [avgCondTableYData;...
                        mean(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                elseif strcmp(whichMetric,'useMedian')
                    errorbar(mean(xCondAvg),median(yCondAvg(~isinf(yCondAvg)),'omitnan'),...
                        std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom,...
                        'o','MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','LineWidth',2);
                    avgCondTableYData = [avgCondTableYData;...
                        median(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                else
                    disp('whichMetric needs to be set to either ''useMean'' or ''useMedian''');
                end
                avgCondTableCondition_Name = [avgCondTableCondition_Name;...
                    obj.conditionNames(find(obj.condIndices == curCond,1))];
                avgCondTableXData = [avgCondTableXData;...
                    mean(xCondAvg,'omitnan')];
                avgCondTableSTDData = [avgCondTableSTDData;...
                    std(yCondAvg(~isinf(yCondAvg)),'omitnan')];
                avgCondTableSEMData = [avgCondTableSEMData;...
                    std(yCondAvg(~isinf(yCondAvg)),'omitnan')/nDenom];
                avgCondTableNData = [avgCondTableNData;...
                    sum(~isnan(yCondAvg) & ~isinf(yCondAvg))];
            else
                avgCondTableCondition_Name = [avgCondTableCondition_Name;...
                    obj.conditionNames(find(obj.condIndices == curCond,1))];
                avgCondTableXData = [avgCondTableXData;NaN];
                avgCondTableYData = [avgCondTableYData;NaN];
                avgCondTableSTDData = [avgCondTableSTDData;NaN];
                avgCondTableSEMData = [avgCondTableSEMData;NaN];
                avgCondTableNData = [avgCondTableNData;0];
            end
                                    
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            if isinf(yMin) && yMin >= 0
                if isinf(yMax) && yMax <= 0
                    yMin = 0;
                    yMax = 1;
                end
            else
                yMin = min([0,1.2*yMin]);
                yMax = 1.2*yMax;
                if yMin == yMax
                    if yMin == 0
                        yMin = yMin - 1;
                        yMax = yMax + 1;
                    else
                        yMin = yMin - 0.1*yMin;
                        yMax = yMax + 0.1*yMin;
                    end
                end
            end
             
            ylim([yMin,yMax]);
            ylabel(strrep(metricName,'_','\_'));
            grid on

            % build full table
            clustTable = table(cTableCondition_Name,cTableSample,cTableFileName,cTableEggChamberID,...
                cTableEggChamberStage,cTableNucleusLabel,cTableCLustLabel, ...
                cTableXData,cTableYData,'VariableNames',{...
                'Condition','sample','sample_InputFileName','eggChamber','eggChamber_Stage',...
                'nuc_Label','clust_Label','x',metricName}); 

            % build table of nuclei averages
            avgNucClustTable = grpstats(clustTable,...
                ["Condition","sample","sample_InputFileName","eggChamber","eggChamber_Stage","nuc_Label"],...
                ["mean","std"],...
                "DataVars",metricName);
            avgNucClustTable = renamevars(avgNucClustTable,{'GroupCount'},{'nClusters'});
            avgNucClustTable.Properties.RowNames = {}; % remove the new row names created by the grpstats function

            % build table of averages by eggChamber
            avgEcClustTable = table(avgEcTableCondition_Name, avgEcTableSample, avgEcTableFileName,...
                avgEcTableEggChamberID, avgEcTableEggChamberStage, avgEcTableXData, avgEcTableYData, ...
                avgEcTableSTDData, avgEcTableSEMData,avgEcTableNData,'VariableNames',...
                {'Condition','sample','sample_InputFileName','eggChamber','eggChamber_Stage',...
                'x',[metricName,'_EggChamberMean'],[metricName,'_EggChamberSTD'],[metricName,'_EggChamberSEM'],...
                'N_clusters'});

            % build table of averages by condition
            avgCondClustTable = table(avgCondTableCondition_Name, avgCondTableXData, avgCondTableYData, ...
                avgCondTableSTDData, avgCondTableSEMData,avgCondTableNData,'VariableNames',...
                {'Condition','x',[metricName,'_ConditionMean'],[metricName,'_ConditionSTD'],[metricName,'_ConditionSEM'],...
                'N_eggChambers'});

        end
        
        %% scatter plot any cluster metric by nucleus
        function fh = plotClusterMetricByNucleus(obj,...
                avgNucClustTable,metricName,idxData,conditionOrder,eggChamberStagesToInclude,figName,alphaVal,varargin)
            % avgNucClustTable is generated by the other scatter plot
            % functions
            
            %collect figure handle if present   
            userDefinedColorMap = [];
            if ~isempty(varargin)
                fh =varargin{1};
                if numel(varargin)>1
                    xOffset = varargin{2};
                    if numel(varargin)>2
                        userDefinedColorMap = varargin{3};
                    end
                end
            else
                fh = [];
                xOffset = 0;
            end

            if isempty(fh)
                fh = figure('Name',figName);
                hold;
                 % generate color maps
                nc = numel(obj.condIndices);
                if isempty(userDefinedColorMap)
                    cmData =cbrewer('qual', 'Dark2', max(nc,3));
                    if nc ==2
                        cmData = cmData([1,3],:);
                    end
                else
                    cmData = userDefinedColorMap;
                end 
            else
                % generate color maps
                nc = numel(obj.condIndices);
                if isempty(userDefinedColorMap)
                    cmData =cbrewer('qual', 'Dark2', max(nc,3));
                    if nc ==2
                        cmData = cmData([1,3],:);
                    end
                else
                    cmData = userDefinedColorMap;
                end    
            end

            % filter out undesired egg chamber stages
            idxStage = ismember(avgNucClustTable.eggChamber_Stage,eggChamberStagesToInclude);
            
            xCtr = xOffset;
            xSampleValsVec = [];
            xSampleIDsVec = {};
            yMin = Inf;
            yMax = -Inf;
            for idxCondOrder=1:numel(conditionOrder)
                i = find(ismember(obj.condIndices,conditionOrder(idxCondOrder)));
                if ~isempty(i)
                    
                    yPlot = avgNucClustTable.(metricName)(idxData ...
                        & ismember(avgNucClustTable.Condition,obj.conditionNames(obj.condIndices(i))) ...
                        & idxStage);
                    if ~isempty(yPlot)
                        nNuclei = numel(yPlot);
                        nucSpacing = obj.cloudWidth*obj.spacingUnit / (nNuclei-1);
                        xCtr = xCtr+obj.condSeparator*obj.spacingUnit;
                        xPlot = xCtr - floor(nNuclei/2)*nucSpacing ...
                                        + (0:(nNuclei-1))*nucSpacing;
                        xMeanCond = mean(xPlot);
                        yMeanCond = mean(yPlot);
                        p = scatter(xPlot,yPlot,'o','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                        alpha(p,alphaVal);
                        p2 = scatter(xMeanCond,yMeanCond,'d','MarkerEdgeColor',cmData(i,:),'MarkerFaceColor',cmData(i,:));
                        alpha(p2,0.9);

                        curID = obj.conditionNames{i};
                        curID = strrep(curID,'_','\_');
                        xSampleValsVec = [xSampleValsVec,xCtr];
                        xSampleIDsVec = [xSampleIDsVec,curID];

                        yMax = max(yMax,max(yPlot(:)));
                        yMin = min(yMin,min(yPlot(:)));
                    end
                end
            end

            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            if isinf(yMin) && yMin >= 0
                if isinf(yMax) && yMax <= 0
                    yMin = 0;
                    yMax = 1;
                end
            else
                yMin = min([0,1.2*yMin]);
                yMax = 1.2*yMax;
                if yMin == yMax
                    if yMin == 0
                        yMin = yMin - 1;
                        yMax = yMax + 1;
                    else
                        yMin = yMin - 0.1*yMin;
                        yMax = yMax + 0.1*yMin;
                    end
                end
            end
             
            ylim([yMin,yMax]);
            ylabel(strrep(metricName,'_','\_'));
            grid on
        end

        %% split a variable name into its constituents
        function [prefix,channel,baseName,suffix] = splitVarName(obj,varName)

            prefix = [];
            channel = [];
            baseName = [];
            suffix = [];

            [varType,baseNameRecognized] = getVarType(obj,varName);
            if isempty(varType)
                return
            end
            
            baseName = baseNameRecognized;
            switch varType
                case 'sample'
                   % do nothing since all values initialized to []
                case 'geometry'
                    k = strfind(varName,baseNameRecognized);
                    p = varName(1:k-1);
                    for i=1:numel(obj.prefixList)
                        if contains(p,obj.prefixList{i})
                            prefix = obj.prefixList{i};
                        end
                    end
                case 'channel'
                    k = strfind(varName,baseNameRecognized);
                    p = varName(1:k-1); % p should be something like nuc_C1
                    for i=1:numel(obj.prefixList)
                        if contains(p,obj.prefixList{i})
                            prefix = obj.prefixList{i};
                            % now that we found the prefix in our string (say 'nuc' in 'nucC1_'), we
                            % collect the channel ID by removing both the
                            % prefix and any '_' and 'C'
                            p = strrep(p,prefix,'');
                            p = strrep(p,'_','');
                            p = strrep(p,'C','');
                            channel = str2double(p);
                        end
                    end

                    q = varName(k+length(baseNameRecognized):end);
                    for i=1:numel(obj.suffixList)
                        if contains(q,obj.suffixList{i})
                            suffix = obj.suffixList{i};
                        end
                    end 
            end
        end

        %% build a variable name based on its constituents
        function varName = buildVarName(obj,prefix,channel,baseName,suffix,chosenVarType)
        % chosenVarType can be either 'channel' or 'geom'

            varName = [];
            varType = [];
            if ismember(baseName,obj.sampleVarList) 
                varType = 'sample';
            end

            % check whether the variable name is ambiguous
            if ismember(baseName,obj.geomVarsBaseNameList) && ismember(baseName,obj.channelVarsBaseNameList)
                if ismember(chosenVarType,{'channel','geom'})
                    varType = chosenVarType;
                else
                    disp(['Cannot build variable names because ambiguous base name ',...
                       chosenVarType,' cannot be resolved into either geometry or channel variable.']);
                    return
                end
            else
                if ismember(baseName,obj.geomVarsBaseNameList)
                    varType = 'geom';
                end

                if ismember(baseName,obj.channelVarsBaseNameList)
                    varType = 'channel';
                end
            end

            if isempty(varType)
                disp(['Cannot build variable, basename ',baseName,' not recognized.']);
                return
            end

            switch varType
                case 'sample'
                    varName = baseName;
                case 'geom'
                    if ~ismember(prefix,obj.prefixList)
                        disp(['could not build variable name, prefix ',prefix,' not recognized']);
                        return
                    else
                        varName = [prefix,'_',baseName];
                    end
                case 'channel'
                    if ~ismember(prefix,obj.prefixList)
                        disp(['could not build variable name, prefix ',prefix,' not recognized']);
                        return
                    else
                        % commented this part out so we can build variable
                        % names before the channels are loaded.
                        %obj.channelList = obj.getChannelList;
                        %if ~ismember(channel,obj.channelList)
                        %    disp(['could not build variable name, channel ',num2str(channel),' not recognized']);
                        %    return
                        %end
    
                        if ~ismember(suffix,obj.suffixList)
                            disp(['could not build variable name, suffix ',suffix,' not recognized']);
                            return
                        end
    
                        varName = [prefix,'_C',num2str(channel),baseName,'_',suffix];
                    end
            end        
        end

end


    methods(Access = 'private')

        %% generate x coordinates by sample/condition for scatter plots
        function [xSampleVals, xSampleIDs] = getSampleXValuesBySample(obj)
            % subfunction used to generate x coordinates for plotting,
            % grouped by sample, then condition.
            % xSampleVals{i}(j) is the x coordinate for condition i, sample j 
            % and matching tick names xSampleIDs{i}{j}
            
            xSampleVals = cell(obj.nConditions,1);
            xSampleIDs = cell(obj.nConditions,1);
            curX = 0;
            for i=1:obj.nConditions
                xSampleVals{i} = zeros(obj.nSamples(i),1);
                xSampleIDs{i} = cell(obj.nSamples(i),1);
                if i>1
                    curX = curX + obj.condSeparator* obj.spacingUnit;
                end
                
                xSampleVals{i}(1:obj.nSamples(i)) = curX + ...
                    (1:obj.nSamples(i)) * obj.sampleSeparator* obj.spacingUnit;

                curX = curX + obj.nSamples(i)*obj.sampleSeparator* obj.spacingUnit;
            
                for j=1:obj.nSamples(i)
                    xSampleIDs{i}{j} = [obj.conditionNames{i},' sample',num2str(j)];
                end
            end
        end
    
        %% generate x coordinates by sample/condition/egg chamber for scatter plots
        function [xEggChamberVals, xEggChamberIDs, xEggChamberStages] = ...
                getSampleXValuesByEggChamberGatedByIdx(obj,eggChamberStagesToInclude,t,idx)
            % subfunction used to generate x coordinates for plotting,
            % grouped by egg chamber, then sample, then condition.
            % xEggChamberVals{i}{j}(k) is the x coordinate for condition i, sample j , egg chamber k
                % with matching tick names xEggChamberIDs{i}{j}{k}
                % and stages xEggChamberStages{i}{j}(k)
            % eggChamberStagesToInclude is a numerical array of stage numbers
                % that should be plotted (e.g. [6,7,8]). 
                % Nuclei/Clusters for other egg chambers will be ignored.
                % You can also set eggChamberStagesToInclude to 'all' as a
                % shortcut to plot all stages.
            % t is whether you want to plot data from the clusters or
                % nuclei table.
            % idx is a logical vector which size matches the number of rows
                % of nucFullT or clustT (depending on what t is set to), which
                % is 1 for the rows to include in the plot, 0 otherwise.
            
            % check that egg chamber data is present
            if isempty(obj.eggChamberNumber) || isempty(obj.eggChamberIDs) || isempty(obj.eggChamberStages)     
                disp(['Variable eggChamberNumber or eggChamberIDs uninitialized',...
                    'or eggChamberStages uninitialized, ',...
                    'cannot generate sample X values by egg chamber']);
                xEggChamberVals = {};
                xEggChamberIDs = {};
                return
            end
            
            % if shorthand 'all' was used for cell stages to include, replace it by
            % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % loop through conditions, samples and egg chambers and
            % generate incremental x values with matching ticks.
            xEggChamberVals = cell(obj.nConditions,1);
            xEggChamberIDs = cell(obj.nConditions,1);
            xEggChamberStages = cell(obj.nConditions,1);
            curX = 0;
            for i=1:obj.nConditions
                xEggChamberVals{i} = cell(obj.nSamples(i),1);
                xEggChamberIDs{i} = cell(obj.nSamples(i),1);
                xEggChamberStages{i} = cell(obj.nSamples(i),1);
                if i>1 && sum((t.condIdx(idx)== obj.condIndices(i)) ~= 0)
                    curX = curX + obj.condSeparator * obj.spacingUnit;
                end
                       
                for j=1:obj.nSamples(i)
                    if j>1 && sum( ...
                            ((t.condIdx(idx)== obj.condIndices(i)) ...
                            & (t.sampleIdx(idx)== obj.sampleIndices{i}(j)) ...
                            )~= 0)
                        curX = curX + obj.sampleSeparator* obj.spacingUnit;    
                    end
                    ctr = 0;

                    xEggChamberVals{i}{j} = zeros(obj.eggChamberNumber{i}(j),1);
                    xEggChamberIDs{i}{j} = cell(obj.eggChamberNumber{i}(j),1);
                    xEggChamberStages{i}{j} = zeros(obj.eggChamberNumber{i}(j),1);
                    
                    for k=1:obj.eggChamberNumber{i}(j)
                        xEggChamberStages{i}{j}(k) = obj.eggChamberStages{i}{j}(k);
                        % only include egg chamber if it isn't part of the
                        % eggChamberStagesToExclude list and there are data
                        % points
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude) ...
                                && sum( ...
                                    ((t.condIdx(idx)== obj.condIndices(i)) ...
                                    & (t.sampleIdx(idx)== obj.sampleIndices{i}(j)) ...
                                    & (t.eggChamberIdx(idx)== obj.eggChamberIDs{i}{j}(k)) ...
                                    )~= 0)
                            if ctr ~=0
                                curX = curX + obj.ecSeparator* obj.spacingUnit;    
                            end
                            xEggChamberVals{i}{j}(k) = curX;
                            xEggChamberIDs{i}{j}{k} = ...
                                [obj.conditionNames{i},' sample',num2str(j),...
                                ' ec',num2str(obj.eggChamberIDs{i}{j}(k))];
                            ctr = ctr+1;
                        else
                            xEggChamberVals{i}{j}(k) = NaN;
                            xEggChamberIDs{i}{j}{k} = [obj.conditionNames{i},' sample',num2str(j),...
                                ' ec',num2str(obj.eggChamberIDs{i}{j}(k))];
                        end
                    end           
                end
            end
        end

        %% generate x coordinates by sample/condition/egg chamber for scatter plots
        function [xEggChamberVals, xEggChamberIDs, xEggChamberStages] = ...
                getSampleXValuesByEggChamber(obj,eggChamberStagesToInclude)
            % subfunction used to generate x coordinates for plotting,
            % grouped by egg chamber, then sample, then condition.
            % xEggChamberVals{i}{j}(k) is the x coordinate for condition i, sample j , egg chamber k
            % with matching tick names xEggChamberIDs{i}{j}{k}
            % and stages xEggChamberStages{i}{j}(k)
            % eggChamberStagesToInclude is a numerical array of stage numbers
            % that should be plotted (e.g. [6,7,8]). 
            % Nuclei for other egg chambers will be ignored.
            % You can also set eggChamberStagesToInclude to 'all' as a
            % shortcut to plot all stages.
            
            % check that egg chamber data is present
            if isempty(obj.eggChamberNumber) || isempty(obj.eggChamberIDs) || isempty(obj.eggChamberStages)     
                disp(['Variable eggChamberNumber or eggChamberIDs uninitialized',...
                    'or eggChamberStages uninitialized, ',...
                    'cannot generate sample X values by egg chamber']);
                xEggChamberVals = {};
                xEggChamberIDs = {};
                return
            end
            
            % if shorthand 'all' was used for cell stages to include, replace it by
            % list of stage numbers.
            if isa(eggChamberStagesToInclude,'char') || isa(eggChamberStagesToInclude,'string')
                if strcmp(eggChamberStagesToInclude,'all')
                    eggChamberStagesToInclude = [0,1,2,3,4,5,6,7,8,9,10];
                end
            end
            
            % loop through conditions, samples and egg chambers and
            % generate incremental x values with matching ticks.
            xEggChamberVals = cell(obj.nConditions,1);
            xEggChamberIDs = cell(obj.nConditions,1);
            xEggChamberStages = cell(obj.nConditions,1);
            curX = 0;
            for i=1:obj.nConditions
                xEggChamberVals{i} = cell(obj.nSamples(i),1);
                xEggChamberIDs{i} = cell(obj.nSamples(i),1);
                xEggChamberStages{i} = cell(obj.nSamples(i),1);
                if i>1
                    curX = curX + obj.condSeparator * obj.spacingUnit;
                end
                       
                for j=1:obj.nSamples(i)
                    if j>1
                        curX = curX + obj.sampleSeparator* obj.spacingUnit;    
                    end
                    ctr = 0;

                    xEggChamberVals{i}{j} = zeros(obj.eggChamberNumber{i}(j),1);
                    xEggChamberIDs{i}{j} = cell(obj.eggChamberNumber{i}(j),1);
                    xEggChamberStages{i}{j} = zeros(obj.eggChamberNumber{i}(j),1);
                    
                    for k=1:obj.eggChamberNumber{i}(j)
                        xEggChamberStages{i}{j}(k) = obj.eggChamberStages{i}{j}(k);
                        % only include egg chamber if it isn't part of the
                        % eggChamberStagesToExclude list
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude)
                            if ctr ~=0
                                curX = curX + obj.ecSeparator* obj.spacingUnit;    
                            end
                            xEggChamberVals{i}{j}(k) = curX;
                            xEggChamberIDs{i}{j}{k} = ...
                                [obj.conditionNames{i},' sample',num2str(j),...
                                ' ec',num2str(obj.eggChamberIDs{i}{j}(k))];
                            ctr = ctr+1;
                        else
                            xEggChamberVals{i}{j}(k) = NaN;
                            xEggChamberIDs{i}{j}{k} = [obj.conditionNames{i},' sample',num2str(j),...
                                ' ec',num2str(obj.eggChamberIDs{i}{j}(k))];
                        end
                    end           
                end
            end
        end

        %% collect names of samples and conditions folders
        function obj = collectConditionsAndSamples(obj,varargin)
            if ~isempty(varargin)
                obj.inFolder = varargin{1};
            end
            disp('Finding out how many conditions and samples there are in the dataset...');
            % collect conditions, i.e. list of subfolders from input folder
            d = dir(obj.inFolder);
            dFolders = d([d(:).isdir]);
            dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
            obj.conditionNames = {dFolders(:).name};
            obj.conditionNames = reshape(obj.conditionNames,...
                numel(obj.conditionNames),1);
            obj.nConditions = numel(dFolders);
            obj.condIndices = 1:obj.nConditions;
            if obj.nConditions > 0
                disp(['found ',num2str(obj.nConditions),' conditions in input folder ',obj.inFolder,'; looking for samples...']);
                % for each condition, collect samples, i.e. list of subfolders
                % within each condition folder:
                for i=1:obj.nConditions
                    d = dir(fullfile(obj.inFolder, obj.conditionNames{i}));
                    dFolders = d([d(:).isdir]);
                    dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
                    obj.sampleNames{i,1} = {dFolders(:).name};
                    obj.sampleNames{i,1} = reshape(obj.sampleNames{i},...
                        numel(obj.sampleNames{i}),1);
                    obj.nSamples(i,1) = numel(dFolders);
                    obj.sampleIndices{i,1} = 1:obj.nSamples(i,1);
                    disp(['    Condition ',num2str(i),'; found ',num2str(obj.nSamples(i,1)),' samples in subfolder ',obj.conditionNames{i},'...']);
                end
            else
                disp(['Did not find any conditions in input folder ',obj.inFolder,'; check that there is data in the input folder.']);
            end
        end
    
        %% get the condition Idx from the name of a condition folder 
        function conditionIdx = getConditionIdx(obj,cName)
            conditionIdx = find(ismember( obj.conditionNames,cName));
        end
        
        %% get the condition and sample Idx from the name of a sample folder 
        % (which should be the name of the initial image)
        function [sample_Idx,conditionIdx] = getSampleIdx(obj,cName,sName)
            conditionIdx = find(ismember( obj.conditionNames,cName));
            sample_Idx = find(ismember( obj.sampleNames{conditionIdx},sName));
        end
    
        %% get the list of expected files to be generated by a complete
        % analysis of the eggChamber data (only channel-independent files).
        function fNames = expectedEggChamberAnalysisFilesChannelIndependentFiles(obj)
            fNames = obj.sampleFOVWiseFileNames(...
                ~obj.sampleFOVWiseChannelDependent);
        end
    
        %% get the list of expected files to be generated by the
        % analysis of the eggChamber data for a specific channel .
        function fNames = expectedEggChamberAnalysisFilesChannelFiles(obj,channelIdx)
           
            nameList = obj.sampleFOVWiseFileNames(...
                obj.sampleFOVWiseChannelDependent);
    
            fNames = {};
            for i=1:numel( nameList )
                fNames = [fNames; ['C',num2str(channelIdx),nameList{i}] ];
            end
        end
        
        %% get the complete list of expected files to be generated by a complete
        % analysis of the eggChamber data. All channels + channel
        % independent files.
        function fNames = expectedEggChamberAnalysisFilesFullList(obj,nChannels)
            fNames = ...
                obj.expectedEggChamberAnalysisFilesChannelIndependentFiles;
    
            for i=1:nChannels
                fNames = [fNames; ...
                    obj.expectedEggChamberAnalysisFilesChannelFiles(i)];
            end
        end
        
        %% get the name of the folder holding the eggChamber related csv files.
        function analysisDir = getEggChamberAnalysisDir(obj,conditionIdx,sample_Idx)
            if isempty(obj.conditionNames)
                disp('condition folder appears empty.');
                analysisDir = [];
                return ;
            end
            analysisDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sample_Idx},...
                obj.sampleFOVCsvFolderName);
        end
        
        %% get the name of the folder holding the eggChamber related segmentation files.
        function eggSegDir = getEggChamberSegDir(obj,conditionIdx,sample_Idx)
            if isempty(obj.conditionNames)
                disp('condition appears empty.');
                eggSegDir = [];
                return ;
            end
            eggSegDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sample_Idx},...
                obj.eggChamberSegFolderName);
        end

        %% get the name of the folder holding the cluster related csv files.
        function clustDir = getClusterDir(obj,conditionIdx,sample_Idx)
            if isempty(obj.conditionNames)
                disp('condition appears empty.');
                clustDir = [];
                return ;
            end
            clustDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sample_Idx},...
                obj.clusterFolderName);
        end
        
        %% check whether all files expected from the analysis of the 
        % eggchamber data are present in the analysis folder
        function [channelsFound,eggChamberSegFound] = areEggChamberAnalysisFilesPresent(obj,...
                                conditionIdx,sample_Idx)
    
            % check that the analysis dir is present
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sample_Idx);
            
            if ~exist(analysisDir, 'dir')
                channelsFound = 0;
                eggChamberSegFound = 0;
                disp(['warning: sample dir ',analysisDir,' not found.']);
                return;
            end
    
            % check that the egg chamber seg dir is present
            eggSegDir = obj.getEggChamberSegDir(conditionIdx,sample_Idx);
            
            if ~exist(eggSegDir, 'dir')
                eggChamberSegFound = 0;
                disp(['optional egg chamber dir ',eggSegDir,' not found.']);
            else
                if exist(fullfile(eggSegDir,obj.eggChamberSegFileName), 'file')
                    eggChamberSegFound = 1;
                else
                    eggChamberSegFound = 0;
                    disp(['optional egg chamber analysis file ',obj.eggChamberSegFileName,' not found.']);
                end
            end
            
            % check that channel-independent analysis files exist
            fNames = ...
                obj.expectedEggChamberAnalysisFilesChannelIndependentFiles;
            flag = 1;
            for i=1:numel(fNames)
                if ~exist(fullfile(analysisDir,fNames{i}),'file')
                       disp(['warning: file ',fullfile(analysisDir,fNames{i}),' not found.']);
                       flag = 0;
                end
            end
            if flag == 0
                channelsFound = 0;
                return;
            end
    
            % check that channel analysis files exist
            channelsFound = 0;
            for i=1:obj.maxChannels
                fNames = ...
                    obj.expectedEggChamberAnalysisFilesChannelFiles(i);
                flag = 1;
                for j=1:numel(fNames)
                    if ~exist(fullfile(analysisDir,fNames{j}),'file')
                           %disp(['warning: file ',fullfile(analysisDir,fNames{j}),' not found.']);
                           flag = 0;
                    end
                end
                if flag == 1
                    channelsFound = i;
                end
            end
        end
        
        %% load the tables corresponding to one sample
        % conditionIDx, sample_Idx relative to conditionNames and
        % sampleNames cell arrays. 
        % loads whole-image metrics (wholeImg and eggChamber prefix),
        % as well as whole nucleus geometry and intensity metrics.
        % the script will load the egg chamber stage info if present,
        % prompting the user for the channel to use as the egg Chamber ID.
        % the last channel should be the correct one based on how the 
        % egg chamber ID is stored by the fiji script.
        
        % outputs a table t that has metrics as variables and individual
        % nuclei as rows.

        % removeNeighbors is a flag - set to 1 to remove variables that contain the word neighbor
        % (these are generated by Fiji but tend to not be useful for our
        % analyses).
        function t = loadEggChamberData(obj,conditionIdx,sample_Idx,removeNeighbors)
            
            % check that data is present
            [c,eggSegFound] = obj.areEggChamberAnalysisFilesPresent(...
                conditionIdx,sample_Idx);
            if c == 0
                t = [];
                disp(['Could not find analysis files for condition ',num2str(obj.condIndices(conditionIdx)),...
                    ' sample ',num2str(obj.sampleIndices{conditionIdx}(sample_Idx))]);
                return
            end
            
            % if there is an egg chamber segmentation file and the channel isnt set yet (i.e. defaulted to zero),
            % prompt the user for the channel that holds the egg Chamber ID,
            % and store it under eggChamberSegChannel property.
            if eggSegFound && (obj.eggChamberSegChannel == 0)
                % assume the egg chamber segmentation channel is the last
                % channel
                prompt = ['Egg Chamber Stage File found - Should I use the last channel ('...
                                ,num2str(c),') as the egg chamber ID? Y/N [Y]: '];
                txt = input(prompt,"s");
                if isempty(txt)
                    txt = 'Y';
                end
                if strcmpi(txt,'Y')
                    obj.eggChamberSegChannel = c;
                else
                    prompt = ['What is the egg chamber ID channel (integer 1-',num2str(c),') ? '];
                    
                    if(obj.eggChamberSegChannel<1 || obj.eggChamberSegChannel>c ...
                            || (floor(obj.eggChamberSegChannel) ~= obj.eggChamberSegChannel))
                        disp('Channel needs to be an integer, using last channel instead');
                        obj.eggChamberSegChannel = c;
                    else
                        obj.eggChamberSegChannel = input(prompt);
                    end
                end
            end
    
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sample_Idx);
            
            tS = []; % full table place holder for data that is single row 
            % - e.g. whole image intensity, sampleROI intensity.
            tE = []; % full table place holder for data that is one row per egg chamber
            % - e.g. whole egg chamber intensity.
            t = []; % full table place holder for data that is mutli row
            % i.e. nucleus-specific 
            % loop through files, load each as a table curT and append the 
            % variables curT holds as new columns into t (or tS if single
            % row)
            for i=1:numel(obj.sampleFOVWiseFileNames)
                if ~obj.sampleFOVWiseChannelDependent(i)
                    % load table
                    curT = readtable( fullfile(analysisDir,...
                        obj.sampleFOVWiseFileNames{i}) );
                    
                    % remove useless variable name (if present)
                    if ismember('Var1',curT.Properties.VariableNames)
                        curT = removevars( curT,'Var1');
                    end
    
                    % add prefix and suffix to variable name, remove
                    % underscores in metric name
                    for k=1:numel(curT.Properties.VariableNames)
                        if strcmp(curT.Properties.VariableNames{k},'InputFileName')
                            curT.Properties.VariableNames{k} = 'sample_InputFileName';

                        elseif ~strcmp(curT.Properties.VariableNames{k},'Label')
                            curT.Properties.VariableNames{k} = ...
                                obj.buildVarName(obj.sampleFOVWisePrefix{i},0,...
                                strrep(curT.Properties.VariableNames{k},'_',''),...
                                obj.eggChamberWiseSuffix{i},'geom');
                        end
                    end
                    
                    % append to full table
                    if obj.sampleFOVWiseSingleRow(i) == 1 % single row table
                        if isempty(tS)
                            tS = curT;
                        else
                            tS = join(tS,curT,'Keys','Label');
                        end
                    elseif obj.sampleFOVWiseSingleRow(i) == 0 % one row per nucleus
                        if isempty(t)
                            t = curT;
                        else
                            t = join(t,curT,'Keys','Label');
                        end
                    elseif obj.sampleFOVWiseSingleRow(i) == 2 % one row per egg chamber
                        if isempty(tE)
                            tE = curT;
                        else
                            tE = join(tE,curT,'Keys','Label');
                        end
                    end
                else
                    for j=1:c
                        % load table
                        curT = readtable( fullfile(analysisDir,...
                            ['C',num2str(j),...
                            obj.sampleFOVWiseFileNames{i}]) );
                        
                        % remove useless variable name (if present)
                        if ismember('Var1',curT.Properties.VariableNames)
                            curT = removevars( curT,'Var1');
                        end
    
                        % add prefix and suffix to variable name
                        for k=1:numel(curT.Properties.VariableNames)
                            if ~strcmp(curT.Properties.VariableNames{k},'Label')
                                curT.Properties.VariableNames{k} = ...
                                    obj.buildVarName(obj.sampleFOVWisePrefix{i},j,...
                                    curT.Properties.VariableNames{k},...
                                    obj.eggChamberWiseSuffix{i},'channel');
                            end
                        end
    
                        % append to full table
                        if obj.sampleFOVWiseSingleRow(i) == 1
                            if isempty(tS)
                                tS = curT;
                            else
                                tS = join(tS,curT,'Keys','Label');
                            end
                        elseif obj.sampleFOVWiseSingleRow(i) == 0
                            if isempty(t)
                                t = curT;
                            else
                                t = join(t,curT,'Keys','Label');
                            end
                        elseif obj.sampleFOVWiseSingleRow(i) == 2
                            if isempty(tE)
                                tE = curT;
                            else
                                tE = join(tE,curT,'Keys','Label');
                            end
                        end
                    end
                end
            end
    
            % combine single row table with nucleus-list table
            tS = repmat(tS,size(t,1),1);
            tS.Label = t.Label;
            t = join(t,tS,'Keys','Label');
    
            % rename Label variable to nuc_Label
            t =renamevars(t,{'Label'},{'nuc_Label'});
    
            % add condition and sample index variables
            cIdx = repmat(conditionIdx,size(t,1),1);
            sIdx = repmat(sample_Idx,size(t,1),1);
            t = addvars(t,cIdx,sIdx,...
                'NewVariableNames',{'cond_Idx','sample_Idx'},...
                'Before',1);
    
            % if egg chamber stage file is present, add a variable for the egg
            % chamber stage.
            if eggSegFound && (obj.eggChamberSegChannel ~= 0)
                
                % get the name of the variable holding the median nucleus
                % intensity in the channel that holds the egg chamber
                % segmenentation ID (the median int will be used as the ID).
                eggChamberIDVariable = obj.buildVarName(...
                    obj.sampleFOVWisePrefix{2},...
                    obj.eggChamberSegChannel,...
                    'Median',...
                    obj.eggChamberWiseSuffix{2},'channel');
                
                % convert the egg chamber ID into and egg chamber stage based
                % on the key saved in the file.
                ecStage = obj.getEggChamberStageData(t.(eggChamberIDVariable),conditionIdx,sample_Idx);
                
                % copy the eggChamber variable with ean easy to interpret file
                % name
                t = addvars(t,t.(eggChamberIDVariable),...
                'NewVariableNames',{'eggChamber_Idx'});
    
                % add egg Chamber stage variable to table
                t = addvars(t,ecStage,...
                'NewVariableNames',{'eggChamber_Stage'});

                % join eggChamber specific metrics to table
                tEout = [];
                isTEoutFilled = 0;
                rowsMissingAtBeginning = 0;
                for i=1:size(t,1)
                    curEC = t.(eggChamberIDVariable)(i); % find egg chamber ID of current nucleus
                    curRow = tE(tE.Label == curEC,:); % select matching row in loaded egg chamber data table
                    
                    if isempty(curRow)
                        curEC = t.(eggChamberIDVariable)(i); % find egg chamber ID of current nucleus
                        
                        if isTEoutFilled == 0
                            rowsMissingAtBeginning = rowsMissingAtBeginning + 1;
                        else
                            for col=1:size(tEout,2)
                                if isa(tEout{1,col},'numeric')
                                    if strcmp(tEout.Properties.VariableNames{col},'Label')
                                        blankRow{1,col} = curEC;
                                    else
                                        blankRow{1,col} = NaN;
                                    end
                                else
                                    blankRow(1,col) = {''};
                                end
                            end
                            tEout = [tEout;blankRow];
                        end
                    else
                        if isTEoutFilled == 0
                            isTEoutFilled = 1;
                            if rowsMissingAtBeginning ==0
                                tEout = curRow;
                            else
                                blankRow = curRow;
                                for col=1:size(curRow,2)
                                    if isa(curRow{1,col},'numeric')
                                        blankRow{1,col} = NaN;
                                    else
                                        blankRow(1,col) = {''};
                                    end
                                end
                                tEout = [repmat(blankRow,rowsMissingAtBeginning,1);curRow];
                            end
                        else
                            tEout = [tEout;curRow];
                        end
                    end
                end
                tEout = addvars(tEout,t.nuc_Label,'NewVariableNames',{'nuc_Label'});
                tEout =renamevars(tEout,{'Label'},{'eggChamber_Label'});
                if ismember('sample_InputFileName',tEout.Properties.VariableNames)
                    tEout = removevars(tEout,{'sample_InputFileName'});
                end
                t = join(t,tEout,'Keys','nuc_Label');
                if ~sum(t.eggChamber_Idx == t.eggChamber_Label) == size(t,1)
                    disp('Error: eggChamber labels loaded are in conflict.');
                end
            end
    
            % remove variables that contain the word neighbor
            if removeNeighbors
                t = obj.removeNeighborVarsFromTable(t);
            end
        end
        
        %% concatenate tables vertically. If variables are missing in one of the
        % tables, it is added with intialized defaulted values to NaN or ''.
        % the order matches that of the variables in t1.
        function tOut = combineEcTables(~,t1, t2)
    
            if isempty(t1)
                tOut = t2;
                return;
            end
    
            if isempty(t2)
                tOut = t1;
                return;
            end
    
            % variable list for each table
            v1 = t1.Properties.VariableNames;
            v2 = t2.Properties.VariableNames;
            
            % list of numeric variables
            s = vartype('numeric');
            s1 = t1(:,s);
            v1num = s1.Properties.VariableNames;
            s2 = t2(:,s);
            v2num = s2.Properties.VariableNames;
            
            % add missing variables to each table so both tables have the same
            % complement of variables, order to match the order of the variables
            % in the other table - the order might break if more than one variable is
            % missing.
            v12 = setdiff(v1,v2);
            if ~isempty(v12)
                for i=1:numel(v12)
                    if ismember(v12{i},v1num)
                        newVar = NaN*zeros(size(t2,1),1);
                    else
                        newVar = repmat({''},size(t2,1),1);
                    end
        
                    idx = find(ismember(v1,v12{i}));
                    if idx == 1
                        if ismember(v1{2},v2)
                            t2 = addvars(t2,newVar,'NewVariableNames',v12(i),'Before',v1{2});
                        else
                            t2 = addvars(t2,newVar,'NewVariableNames',v12(i));
                        end
                        v2 = t2.Properties.VariableNames;
                    else
                        if ismember(v1{idx-1},v2)
                            t2 = addvars(t2,newVar,'NewVariableNames',v12(i),'After',v1{idx-1});
                        else
                            t2 = addvars(t2,newVar,'NewVariableNames',v12(i));
                        end
                        v2 = t2.Properties.VariableNames;
                    end          
                end
            end
        
            v21 = setdiff(v2,v1);
            if ~isempty(v21)
                for i=1:numel(v21)
                    if ismember(v21{i},v2num)
                        newVar = NaN*zeros(size(t1,1),1);
                    else
                        newVar = repmat({''},size(t1,1),1);
                    end
        
                    idx = find(ismember(v2,v21{i}));
                    if idx == 1
                        if ismember(v2{2},v1)
                            t1 = addvars(t1,newVar,'NewVariableNames',v21(i),'Before',v2{2});
                        else
                            t1 = addvars(t1,newVar,'NewVariableNames',v21(i));
                        end
                        v1 = t1.Properties.VariableNames;
                    else
                        if ismember(v2{idx-1},v1)
                            t1 = addvars(t1,newVar,'NewVariableNames',v21(i),'After',v2{idx-1});
                        else
                            t1 = addvars(t1,newVar,'NewVariableNames',v21(i));
                        end
                        v1 = t1.Properties.VariableNames;
                    end          
                end
            end
        
            % re-order variables in t2 to match the order in t1 if needed
            v1 = t1.Properties.VariableNames;
            v2 = t2.Properties.VariableNames;
        
            idx12 = zeros(numel(v1,1));
            for i=1:numel(v1)
                idx12(i) = find(ismember(v2,v1(i)));
            end
        
            t2 = t2(:,idx12);
            tOut = [t1;t2]; 
        end
        
        %%
        function ecStage = getEggChamberStageData(obj, eggChamber_Idx,conditionIdx,sample_Idx)
            % load cell stage key file
            eggSegDir = obj.getEggChamberSegDir(conditionIdx,sample_Idx);
            ecStageFile = fullfile(eggSegDir,obj.eggChamberSegFileName);
    
            % reads the egg chamber stage key file as a 2-column table w
            % variable names eggChamber_Idx and eggChamber_Stage.
            tEC = readtable(ecStageFile);
            tEC = renamevars(tEC,{'eggChamberID','eggChamberStage'},{'eggChamber_Idx','eggChamber_Stage'});

            % convert key file to Containers.map object.
            ecMap = containers.Map(tEC.eggChamber_Idx,tEC.eggChamber_Stage);
            ecStage = zeros(numel(eggChamber_Idx),1);
            ecVals = unique(eggChamber_Idx);
            for i=1:numel(ecVals)
                if ~isKey(ecMap,ecVals(i))
                    ecStage(eggChamber_Idx==ecVals(i)) = 0;
                else
                    ecStage(eggChamber_Idx==ecVals(i)) = cell2mat(values(ecMap,{ecVals(i)}));
                end
            end
        end

        %%
        function reOrderNucleiVariables(obj)
            % order should go:
            % input file name
            % cond_Idx
            % sample_Idx
            % eggChamber_Idx (optional)
            % eggchamber stage (optional)
            % nuc_Label
            % nuc_Volume 
            
            vStart = obj.nucT.Properties.VariableNames{1};
            obj.nucT = movevars(obj.nucT,'sample_InputFileName','Before',vStart);
            obj.nucT = movevars(obj.nucT,'cond_Idx','After','sample_InputFileName');
            obj.nucT = movevars(obj.nucT,'sample_Idx','After','cond_Idx');
            if ismember('eggChamber_Idx',obj.nucT.Properties.VariableNames)
                obj.nucT = movevars(obj.nucT,'eggChamber_Idx','After','sample_Idx');
            end

            if ismember('eggChamber_Stage',obj.nucT.Properties.VariableNames)
                obj.nucT = movevars(obj.nucT,'eggChamber_Stage','After','eggChamber_Idx');
                obj.nucT = movevars(obj.nucT,'nuc_Label','After','eggChamber_Stage');
            else
                if ismember('eggChamber_Idx',obj.nucT.Properties.VariableNames)
                    obj.nucT = movevars(obj.nucT,'nuc_Label','After','eggChamber_Idx');
                else
                    obj.nucT = movevars(obj.nucT,'nuc_Label','After','sample_Idx');
                end
            end
            obj.nucT = movevars(obj.nucT,'nuc_Volume','After','nuc_Label');
        end

        %%
        function sortNucleiRowsByEggChamber(obj)
            if ismember('eggChamber_Idx',obj.nucT.Properties.VariableNames)
                obj.nucT = sortrows(obj.nucT,{'cond_Idx','sample_Idx','eggChamber_Idx','nuc_Label'});
            else
                disp('Cannot sort table by egg chamber, eggChamber_Idx column missing. Skipping.');
            end

        end

        %% find the number of channels in the table
        function [channelIndices,nChannels] = getChannelList(obj)
            allVars = obj.nucT.Properties.VariableNames;
            channelIndices = [];
            for i=1:numel(allVars)
                v = obj.getVarType(allVars{i});
                if strcmp(v,'channel')
                    [~,c,~,~] = obj.splitVarName(allVars{i});
                    channelIndices = [channelIndices,c];
                end
            end
            channelIndices = channelIndices(~isnan(channelIndices));
            if obj.eggChamberSegChannel ~=0
                channelIndices = ...
                    union(channelIndices,obj.eggChamberSegChannel);
            end
            channelIndices = unique(channelIndices);
            nChannels = numel(channelIndices);

            obj.channelList = channelIndices;
        end

        

        
        
        %% find whether a given variable (varName) is a sample variable, a
        % geometry variable or a channel variable. Also outputs the
        % baseName of the variable. Returns empty values if unrecognized.
        function [varType,baseNameRecognized]= getVarType(obj,varName)
            if ismember(varName,obj.sampleVarList)
                varType = 'sample';
                baseNameRecognized = ...
                    obj.sampleVarList(ismember(obj.sampleVarList,varName));
                return
            end
                
            % need to test for geometry variables first, otherwise
            % possible ambiguity between MeanBreadth (geometry) and Mean (intensity)
            for i=1:numel(obj.geomVarsBaseNameList)
                if contains(varName,obj.geomVarsBaseNameList{i})
                    varType = 'geometry'; 
                    baseNameRecognized = obj.geomVarsBaseNameList{i};
                    return
                end
            end

            for i=1:numel(obj.channelVarsBaseNameList)
                if contains(varName,obj.channelVarsBaseNameList{i})
                    varType = 'channel'; 
                    baseNameRecognized = obj.channelVarsBaseNameList{i};
                    return
                end
            end
            
            varType = '';
            baseNameRecognized = '';
            disp(['type of variable ',varName,' was not recognized.']);
        end
     
        %% get the IDs and stages of all egg chambers in the experiment.
        % output:
        % ecN is a nested cell array where ecN{i}{j} is the number of egg
            % chambers segmented in condition i sample j.
        % ecID is a nested cell array where ecID{i}{j} is a ecN{i}{j} x 1
            % numerical array holding the IDs of the segmented egg chambers in condition i sample j.
        % ecStage is a nested cell array where ecStage{i}{j} is a ecN{i}{j} x 1
            % numerical array holding the stages of the segmented egg chambers in condition i sample j.
            % defaulted to zero if the egg chamber stage info is missing.
        % ecNumNucPerEC is a nested cell array where ecNumNucPerEC{i}{j} is a ecN{i}{j} x 1
            % numerical array holding the number of nuclei per segmented egg chambers in condition i sample j.
        % ecN, ecID and ecStage are defaulted to {} if the eggChamber_Idx
            % variable is absent from the data table.
        function  [ecID,ecN,ecStage,ecNumNucPerEC] = getEggChamberIDs(obj)
            % check whether the eggChamber_Idx variable is present. If not, return
            % empty cell arrays.
            if ~ismember('eggChamber_Idx',obj.nucT.Properties.VariableNames)
                disp('Cannot compile egg chamber IDs, column missing. Skipping.');
                ecID = {};
                ecN = {};
                ecStage = {};
                ecNumNucPerEC = {};
                return
            end

            % loop through the conditions and samples and collect the IDs
            % and stages for each condition/sample.
            ecID = cell(numel(obj.condIndices),1);
            ecStage = cell(numel(obj.condIndices),1);
            ecNumNucPerEC = cell(numel(obj.condIndices),1);
            ecN = cell(numel(obj.condIndices),1);
            for i=1:numel(obj.condIndices)
                ecID{i} = cell(obj.nSamples(i),1);
                ecN{i} = zeros(obj.nSamples(i),1);
                ecNumNucPerEC{i} = cell(obj.nSamples(i),1);
                ecStage{i} = cell(obj.nSamples(i),1);
                for j=1:obj.nSamples(i)
                    ecID{i}{j} = ...
                        unique(obj.nucT.eggChamber_Idx(...
                        obj.nucT.cond_Idx == i & obj.nucT.sample_Idx == j));
                    ecN{i}(j) = numel(ecID{i}{j});

                    curT = obj.nucT( obj.nucT.cond_Idx == i & obj.nucT.sample_Idx == j,:);
                    for k=1:numel(ecID{i}{j})
                        ecNumNucPerEC{i}{j}(k) = sum(curT.eggChamber_Idx == ecID{i}{j}(k));
                        if ismember('eggChamber_Stage',obj.nucT.Properties.VariableNames)
                            ecStage{i}{j}(k) = ...
                            unique(...
                            curT.eggChamber_Stage(curT.eggChamber_Idx == ecID{i}{j}(k)));
                        else
                            ecStage{i}{j} = zeros(size(ecID));
                        end
                    end
                end
            end
            obj.eggChamberIDs = ecID;
            obj.eggChamberStages = ecStage;
            obj.eggChamberNumber = ecN;
            obj.eggChamberNumNucPerEC = ecNumNucPerEC;
        end      

    %% add nucleus and egg chamber ID to the left of the cluster table
        function clustT = appendNucVarsToClustTable(obj,clustT,cond_Idx,sample_Idx,nuc_Label)
            % nucleus variables to be appended to the left of the cluster
            % table to ease linking clusters to nuclei.

            n = size(clustT,1);
            % looping backwards because we are adding variables to the left of the table
            % this way at the end the order of the added variables in the
            % final table reflects the order in nucVarsToIncludeInClustTable
            for i=numel(obj.nucVarsToIncludeInClustTable):-1:1 
                if ismember(obj.nucVarsToIncludeInClustTable{i},obj.nucT.Properties.VariableNames)
                    curT = repmat(...
                        obj.nucT.(obj.nucVarsToIncludeInClustTable{i})(...
                        obj.nucT.cond_Idx == cond_Idx ...
                        & obj.nucT.sample_Idx == sample_Idx ...
                        & obj.nucT.nuc_Label == nuc_Label),...
                        n,1);
                    clustT = addvars(clustT,curT,...
                        'NewVariableNames',obj.nucVarsToIncludeInClustTable(i),...
                        'Before',clustT.Properties.VariableNames{1});
                end
            end
        end

        %% remove variables containing the word neighbors in them
        function tOut = removeNeighborVarsFromTable(~,tIn)
            idx = cellfun(@contains,tIn.Properties.VariableNames,...
                    repmat( {'Neighbors'},size(tIn.Properties.VariableNames) ),...
                    'UniformOutput',0);
            idx = cell2mat(idx);
            tOut = tIn(:,~idx);
        end

        %% load cluster data for a given nucleus
        % output are two tables, one which is nucleus metrics (nucleoplasm and nucleoli metrics)
        % one which is cluster metrics.
        function [curNucT,curClustT] = loadSampleClustTables(obj,...
                cond_Idx,sample_Idx,nuc_Label,removeNeighbors)

            % where the data is
            clustDir = obj.getClusterDir(cond_Idx,sample_Idx);
            
            curNucT = table();
            curClustT = table();
            for i=1:numel(obj.clusterWiseFileNames)
                if ~obj.clusterWiseChannelDependent(i)
                    curFileName = fullfile(clustDir,...
                        ['nuc',num2str(nuc_Label),obj.clusterWiseFileNames{i}]);
                    % load table
                    if exist(curFileName,"file")
                        curT = readtable(curFileName);
                    else
                        disp(['Could not find file ',curFileName,...
                            ' data will be missing ...']);
                        curT = table();
                    end

                    % remove variables that contain the word neighbor
                    if removeNeighbors
                        curT = obj.removeNeighborVarsFromTable(curT);
                    end

                    % add prefix and suffix to variable name, remove
                    % underscores in the metric name
                    for k=1:numel(curT.Properties.VariableNames)
                        if ~strcmp(curT.Properties.VariableNames{k},'Label')
                            curT.Properties.VariableNames{k} = ...
                                obj.buildVarName(obj.clusterWisePrefix{i},0,...
                                strrep(curT.Properties.VariableNames{k},'_',''),...
                                obj.clusterWiseSuffix{i},'geom');
                        end
                    end

                    tic;
                    % make sure the label of the nucleus table is the
                    % nucleus ID (it can be 255 in some cases)
                    if obj.clusterWiseSingleRow(i)
                        curT.Label = nuc_Label;
                    end

                    % join loaded table to sample-wise table
                    if ~isempty(curT)
                        if obj.clusterWiseSingleRow(i)
                            if isempty(curNucT)
                                curNucT = curT;
                            else
                                curNucT = join(curNucT,curT,'Keys','Label');
                            end
                        else
                            if isempty(curClustT)
                                curClustT = curT;
                            else
                                curClustT = join(curClustT,curT,'Keys','Label');
                            end
                        end
                    end

                else
                    
                    for j=1:numel(obj.channelList)
                        curFileName = fullfile(clustDir,...
                            ['C',num2str(obj.channelList(j)),'_nuc',num2str(nuc_Label),...
                            obj.clusterWiseFileNames{i}]);
                        %disp(curFileName);
                        % load table
                        if exist(curFileName,"file")
                            curT = readtable(curFileName);
                        else
                            disp(['Could not find file ',curFileName,...
                                ' data will be missing ...']);
                            curT = table();
                        end

                        % remove variables that contain the word neighbor
                        if removeNeighbors
                            curT = obj.removeNeighborVarsFromTable(curT);
                        end
                        
                        % add prefix and suffix to each variable name
                        for k=1:numel(curT.Properties.VariableNames)
                            if ~strcmp(curT.Properties.VariableNames{k},'Label')
                                oldName = curT.Properties.VariableNames{k};
                                newName = obj.buildVarName(obj.clusterWisePrefix{i},...
                                    obj.channelList(j),...
                                    strrep(oldName,'_',''),...
                                    obj.clusterWiseSuffix{i},'channel');
                                curT.Properties.VariableNames{k} = newName;
                            end
                        end

                        % make sure the label of the nucleoplasm and nucleoli tables is the
                        % nucleus ID (it can be 255 in some cases)
                        if obj.clusterWiseSingleRow(i)
                            curT.Label = nuc_Label;
                        end

                        % join loaded table to sample-wise table
                        if ~isempty(curT)
                            if obj.clusterWiseSingleRow(i)
                                if isempty(curNucT)
                                    curNucT = curT;
                                else
                                    curNucT = join(curNucT,curT,'Keys','Label');
                                end
                            else
                                if isempty(curClustT)
                                    curClustT = curT;
                                else
                                    curClustT = join(curClustT,curT,'Keys','Label');
                                end
                            end
                        end
                    end
                end    
            end 

            % rename Label variable to nuc_Label/clust_Label
            if ~isempty(curNucT)
                curNucT =renamevars(curNucT,{'Label'},{'nuc_Label'});
            end
            if ~isempty(curClustT)
                curClustT =renamevars(curClustT,{'Label'},{'clust_Label'});
            end

            % in cluster table, copy one of the channel volume variables to
            % a Volume variable 
            % (this avoids all the volume info getting erased when streamlining the table)
            idx = cellfun(@contains,curClustT.Properties.VariableNames,...
                    repmat( {'Volume'},size(curClustT.Properties.VariableNames) ),...
                    'UniformOutput',0);
            idx = cell2mat(idx);
            if sum(idx)>=1
                idx = find(idx,1);
                curClustT = addvars(curClustT,...
                    curClustT.(curClustT.Properties.VariableNames{idx}),...
                    'NewVariableNames',{'clust_Volume'});
            end
        end
        
        %% using nuclei table, collect list of nuclei labels expected in current sample
        function nucIDList = getExpectedNucIDsFromTable(obj,cond_Idx, sample_Idx)
            nucIDList = obj.nucT.nuc_Label(...
                obj.nucT.cond_Idx == cond_Idx ...
                & obj.nucT.sample_Idx == sample_Idx);
        end
    end
end
