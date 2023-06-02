
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
    end

    properties (GetAccess = 'private', SetAccess = 'private')

        % ******************** eggchamber data I/O properties

        % name of the subfolder where the FOV-level (nucleus-wise, sampleROI-wise and wholeImg-wise) data is stored
        sampleFOVCsvFolderName = 'eggChamberCSV';
        % The folder should contain following files:
        % - allNucGeom.csv: geometry metrics of each nucleus
        % - C1_allNucInt.csv: nuclear-wise intensity metrics in channel 1
        % - C1_eggChamberInt.csv: sampleROI-wise intensity metrics in channel 1 (sample ROI is one giant ROI that covers broadly all egg chambers of the image)
        % - C1_wholeImgInt.csv: image-wise intensity metrics in channel 1
        sampleFOVWiseFileNames = {'allNucGeom.csv';...
            '_allNucInt.csv';...
            '_wholeImgInt.csv';...
            '_eggChamberInt.csv'};

        % whether each file name of the previous list is channel-dependent or not
        sampleFOVWiseChannelDependent = logical([0;...
            1;...
            1;...
            1;]);
    
        % whether each file name of the previous list contains a single row (whole image data) or
        % one row per nucleus
        sampleFOVWiseSingleRow = logical([0;...
            0;...
            1;...
            1;]);
        
        % prefix to add in front of each variable name when loading the
        % corresponding table
        sampleFOVWisePrefix ={'nuc';...
            'nuc';...
            'wholeImg';...
            'sampleROI'};

        % name of the optional subfolder where the 2D segmentations of the egg chambers are stored
        eggChamberSegFolderName = 'eggChamberSEG';
        % if used, the folder should hold a file called:
        eggChamberSegFileName = 'eggChamberStages.csv';
    
        % suffix to add at the end of each variable name when loading the
        % corresponding table
        eggChamberWiseSuffix ={'';...
            'raw';...
            'raw';...
            'raw'};
        
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
        'NumClusters'}; 
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
        sampleROISubtractedSuffix = 'sampleROISub';
        wholeImgSubtractedSuffix = 'wholeImgSub';
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
        backgroundIntensityPrefixList;

        % when performing background subtraction, use the median intensity
        % of the reference region.
        metricToSubtract = 'Median';

        % nucleus variables to include in cluster table
        nucVarsToIncludeInClustTable = {'sample_InputFileName','cond_Idx',...
                'sample_Idx','eggChamber_Idx','eggChamber_Stage','nuc_Label'};

        % ************************** plotting settings    

        % used to set the spacing between cloud of points belonging to different nuclei    
        spacingUnit = 1;

        % how much separation there is between conditions in units of spacingUnit
        condSeparator = 3; 

        % how much separation there is between samples in units of spacingUnit
        sampleSeparator = 2;

        % how much separation there is between egg chambers in units of spacingUnit
        ecSeparator = 1.5
        
        % free space between the spot clouds belonging to neighboring
        % nuclei in units of spacingUnit
        freeSpaceBetweenSamples = 0.2;

    end

    methods (Access = 'public')

        %% initialize object from input root folder
        function obj = eggChamberDataSet(inputFolder)
            obj.inFolder = inputFolder;
            obj.collectConditionsAndSamples;

            obj.prefixList = {obj.nucPrefix,obj.wholeImgPrefix,obj.sampleROIPrefix,obj.clusterPrefix,obj.nucleoliPrefix,obj.plasmPrefix};
            obj.suffixList = {obj.rawSuffix,obj.plasmCorrSuffix,obj.sampleROISubtractedSuffix,obj.wholeImgSubtractedSuffix};
            obj.backgroundIntensityPrefixList = {obj.wholeImgPrefix,obj.sampleROIPrefix};

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
            nucStatsVars = setdiff(obj.nucT.Properties.VariableNames, ...
                obj.clustT.Properties.VariableNames);
            nucStatsVars = [nucStatsVars,'nuc_Label'];
            idx = find(ismember(obj.nucT.Properties.VariableNames,nucStatsVars));
            
            tJoin = table();
            for ctr=1:size(obj.clustT,1)
                i = obj.clustT.cond_Idx(ctr);
                j = obj.clustT.sample_Idx(ctr);
                k = obj.clustT.nuc_Label(ctr);

                % extract nuc stats for current nucleus
                nT = obj.nucT(...
                   obj.nucT.cond_Idx == i...
                   & obj.nucT.sample_Idx == j ...
                   & obj.nucT.nuc_Label == k,idx);

                nT2 = join(obj.clustT(ctr,:),nT,'Keys','nuc_Label');
                tJoin = obj.combineEcTables(tJoin,nT2);  
            end
        end

        %% add average cluster stats to nuc tables
        function tJoin = addAverageClusterStatsToNucTable(obj)
            % generate list of cluster variables to average out
            cv = obj.clustT.Properties.VariableNames;
            idx = find(cell2mat( cellfun( @contains, cv,...
                        repmat({'clust'}, size(cv)),...
                        'UniformOutput',0) )); 
            cIdx = find(ismember(cv,{'clust_Label'}));
            idx = setdiff(idx,cIdx); % indices of the variables to average (in clustT table)            
            
            t = table();
            for ctr=1:size(obj.nucT,1)
                i = obj.nucT.cond_Idx(ctr);
                j = obj.nucT.sample_Idx(ctr);
                k = obj.nucT.nuc_Label(ctr);
                
                % extract cluster values for current nucleus
                cT = obj.clustT(...
                   obj.clustT.cond_Idx == i...
                   & obj.clustT.sample_Idx == j ...
                   & obj.clustT.nuc_Label == k,idx);

                if ~isempty(cT)
                    nC = size(cT,1); % number of clusters in current nucleus

                    % average/std cluster metrics and keep only meaningful variables
                    cT1 = varfun(@mean,cT,'InputVariables',@isnumeric);
                    cT2 = varfun(@std,cT,'InputVariables',@isnumeric);
                    
                    % update variable names from mean_clust to avgClust
                    newVars = cellfun(@strrep,cT1.Properties.VariableNames,...
                        repmat({'mean_clust'},size(cT1.Properties.VariableNames)),...
                        repmat({'nucAvgClust'},size(cT1.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT1 = renamevars(cT1,cT1.Properties.VariableNames,newVars);
                    
                    % update variable names from std_clust to stdClust
                    newVars = cellfun(@strrep,cT2.Properties.VariableNames,...
                        repmat({'std_clust'},size(cT2.Properties.VariableNames)),...
                        repmat({'nucStdClust'},size(cT2.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT2 = renamevars(cT2,cT2.Properties.VariableNames,newVars);
                    
                    % add a nuc_Label variable to join tables (and
                    % validation later)
                    cT1 = addvars(cT1,k,'NewVariableNames',{'nuc_Label'});
                    cT2 = addvars(cT2,k,'NewVariableNames',{'nuc_Label'});
                    cT = join(cT1,cT2,'Keys','nuc_Label');

                    % compute the number of clusters per nucleus 
                    cT = addvars(cT,nC,'NewVariableNames',{'nuc_NumClusters'});   
                else
                    % no clusters found in current nucleus, generating a
                    % row of NaN values
                    cT = obj.clustT(1,idx);
                    for cc = 1:size(cT,2)
                        cT{1,cc} = NaN;
                    end

                    % update variable names from clust to avgClust
                    newVars = cellfun(@strrep,cT.Properties.VariableNames,...
                        repmat({'clust'},size(cT.Properties.VariableNames)),...
                        repmat({'avgClust'},size(cT.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT1 = renamevars(cT,cT.Properties.VariableNames,newVars);

                    % same for stdClust
                    newVars = cellfun(@strrep,cT.Properties.VariableNames,...
                        repmat({'clust'},size(cT.Properties.VariableNames)),...
                        repmat({'stdClust'},size(cT.Properties.VariableNames)),...
                        'UniformOutput',0);
                    cT2 = renamevars(cT,cT.Properties.VariableNames,newVars);

                    % add a nuc_Label variable (used for joinig and validation later)
                    cT1 = addvars(cT1,k,'NewVariableNames',{'nuc_Label'});
                    cT2 = addvars(cT2,k,'NewVariableNames',{'nuc_Label'});
                    cT = join(cT1,cT2,'Keys','nuc_Label');

                    % compute the number of clusters per nucleus 
                    cT = addvars(cT,0,'NewVariableNames',{'nuc_NumClusters'});

                end

                % append current nuclei cluters metrics to global table
                t = obj.combineEcTables(t,cT);  
            end
            
            % join t (average cluster metrics per nucleus) with the nucT table
            % holding all other nuclei metrics
            if size(t,1) ~= size(obj.nucT,1)
                disp(['Cannot append average cluster values to nucleus table; ',...
                    'Nuc table has ',num2str(size(obj.nucT,1)),' rows while ',...
                    'Clust table has ',num2str(size(t,1)),' rows.']);
            else
                % add dummy column to each variable to use as 
                dummyKey = (1:size(obj.nucT,1))';
                nuc2 = addvars(obj.nucT,dummyKey,'NewVariableNames',{'Key'});
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
            
            obj.nucT = bgCorrTable(obj.nucT);
            obj.nucFullT = bgCorrTable(obj.nucFullT);

            function tOut = bgCorrTable(tIn)
                tOut = tIn;

                [c,nChannels] = obj.getChannelList;
                
                % collect the background columns to subtract
                varToSubtract = cell(...
                    nChannels,numel(obj.backgroundIntensityPrefixList));
    
                for i=1:nChannels
                    for j=1:numel(obj.backgroundIntensityPrefixList)
                        
                        varToSubtract{i,j} = obj.buildVarName(...
                            obj.backgroundIntensityPrefixList{j},c(i),...
                        obj.metricToSubtract,obj.rawSuffix,'channel');
                        
                    end
                end
                
                % keep only in t the variables that need subtracting
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
                for i=1:numel(obj.backgroundIntensityPrefixList)
                    curIdx = find(cell2mat( cellfun( @contains, varList,...
                        repmat(obj.backgroundIntensityPrefixList(i), size(varList)),...
                        'UniformOutput',0) ));
                    idx = setdiff(idx,curIdx);
                end
                
                varList = tIn.Properties.VariableNames(idx);
                
                % loop through channels and background correct
                for i=1:nChannels
                
                    % find variables in t that are relative to the current color channel
                    idxC = cell2mat( cellfun( @contains, varList,...
                        repmat({['_C',num2str(c(i))]}, size(varList)),...
                        'UniformOutput',0) );
                
                    curVarList = varList(idxC);
                    for j=1:numel(obj.backgroundIntensityPrefixList)
                        for k =1:numel(curVarList)
                            newVarName = strrep(curVarList{k},...
                                ['_',obj.rawSuffix],...
                                ['_',obj.backgroundIntensityPrefixList{j},'Subtracted']);
    
                            newVar = tIn.(curVarList{k}) - tIn.(varToSubtract{i,j});
    
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
                    x = obj.nucT.(varName)(...
                        obj.nucT.cond_Idx ==obj.condIndices(j) ...
                        & obj.nucT.sample_Idx == s(k));
            
                    % generate slightly offset x coordinates for each nucleus,
                    % centered around the sample X
                    nNuclei = size(x,1);
                    if nNuclei >1
                        % spacing between nuclei
                        nucSpacing = obj.spacingUnit ...
                            * (1-2*obj.freeSpaceBetweenSamples)/(nNuclei-1);
            
                        % x coordinate for each nucleus of current condition/sample
                        curXPlot = xSampleVals{j}(k) - floor(nNuclei/2)*nucSpacing ...
                            + (0:(nNuclei-1))*nucSpacing;
            
                        % y coordinate for each nucleus of current condition/sample
                        curYPlot = obj.nucT.(varName)(...
                            obj.nucT.cond_Idx ==obj.condIndices(j) ...
                            & obj.nucT.sample_Idx == s(k))';
            
                    elseif nNuclei == 1
                        curXPlot = xSampleVals{j}(k);
                        curYPlot = obj.nucT.(varName)(...
                            obj.nucT.cond_Idx ==obj.condIndices(j) ...
                            & obj.nucT.sample_Idx == s(k))';
            
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
            ylabel(baseName);

        end

        %% scatter plot a metric by egg chamber
        % prefix: any of 'nuc', 
        function scatterPlotNucleiMetricByEggChamber(obj,prefix,channel,baseName,suffix,eggChamberStagesToInclude)
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
            if ~ismember( varName, obj.nucT.Properties.VariableNames)
                disp(['Variable ',varName,' absent from table, cannot plot.']);
                return
            end

            % build figure
            figure('Name',strrep(varName,'_',' '));
            hold;
            
            % collect the x values to plot each condition/sample at
            [xEggChamberVals, xEggChamberIDs] = ...
                obj.getSampleXValuesByEggChamber(eggChamberStagesToInclude);

            % collect the values of the metric for each condition/sample
            xPlot = [];
            yPlot = [];
            xSampleValsVec = [];
            xSampleIDsVec = {};
            for i=1:numel(obj.condIndices)
                s = obj.sampleIndices{i};
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        if ismember(obj.eggChamberStages{i}{j}(k),eggChamberStagesToInclude)
                            % collect values for the desired metric from all nuclei for the
                            % current sample/condition
                            x = obj.nucT.(varName)(...
                                obj.nucT.cond_Idx ==obj.condIndices(i) ...
                                & obj.nucT.sample_Idx == s(j) ...
                                & obj.nucT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k));
                    
                            % generate slightly offset x coordinates for each nucleus,
                            % centered around the sample X
                            nNuclei = size(x,1);
                            if nNuclei >1
                                % spacing between nuclei
                                nucSpacing = obj.spacingUnit ...
                                    * (1-2*obj.freeSpaceBetweenSamples) / (nNuclei-1);
                    
                                % x coordinate for each nucleus of current condition/sample
                                curXPlot = xEggChamberVals{i}{j}(k) ...
                                    - floor(nNuclei/2)*nucSpacing ...
                                    + (0:(nNuclei-1))*nucSpacing;
                    
                                % y coordinate for each nucleus of current condition/sample
                                curYPlot = obj.nucT.(varName)(...
                                    obj.nucT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucT.sample_Idx == s(j)...
                                    & obj.nucT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 1
                                curXPlot = xEggChamberVals{i}{j}(k);
                                curYPlot = obj.nucT.(varName)(...
                                    obj.nucT.cond_Idx ==obj.condIndices(i) ...
                                    & obj.nucT.sample_Idx == s(j)...
                                    & obj.nucT.eggChamber_Idx == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 0
                                curXPlot = [];
                                curYPlot = [];
                            end
                    
                            % append coordinates of current condition/sample to global list
                            xPlot = [xPlot,curXPlot];
                            yPlot = [yPlot,curYPlot];
                    
                            xSampleValsVec = [xSampleValsVec,xEggChamberVals{i}{j}(k)];
                            xSampleIDsVec = [xSampleIDsVec,xEggChamberIDs{i}{j}{k}];
                        end
                    end
                end
            end
            plot(xPlot,yPlot,'o');
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylabel(baseName);
            grid on
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
    
            % collect conditions, i.e. list of subfolders from input folder
            d = dir(obj.inFolder);
            dFolders = d([d(:).isdir]);
            dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
            obj.conditionNames = {dFolders(:).name};
            obj.conditionNames = reshape(obj.conditionNames,...
                numel(obj.conditionNames),1);
            obj.nConditions = numel(dFolders);
            obj.condIndices = 1:obj.nConditions;

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
                return;
            end
    
            % check that the egg chamber seg dir is present
            eggSegDir = obj.getEggChamberSegDir(conditionIdx,sample_Idx);
            
            if ~exist(eggSegDir, 'dir')
                eggChamberSegFound = 0;
            else
                if exist(fullfile(eggSegDir,obj.eggChamberSegFileName), 'file')
                    eggChamberSegFound = 1;
                else
                    eggChamberSegFound = 0;
                end
            end
            
            % check that channel-independent analysis files exist
            fNames = ...
                obj.expectedEggChamberAnalysisFilesChannelIndependentFiles;
            flag = 1;
            for i=1:numel(fNames)
                if ~exist(fullfile(analysisDir,fNames{i}),'file')
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
            % - e.g. whole image intensity, whole egg chamber intensity.
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
                    if obj.sampleFOVWiseSingleRow(i)
                        if isempty(tS)
                            tS = curT;
                        else
                            tS = join(tS,curT,'Keys','Label');
                        end
                    else
                        if isempty(t)
                            t = curT;
                        else
                            t = join(t,curT,'Keys','Label');
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
                        if obj.sampleFOVWiseSingleRow(i)
                            if isempty(tS)
                                tS = curT;
                            else
                                tS = join(tS,curT,'Keys','Label');
                            end
                        else
                            if isempty(t)
                                t = curT;
                            else
                                t = join(t,curT,'Keys','Label');
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
