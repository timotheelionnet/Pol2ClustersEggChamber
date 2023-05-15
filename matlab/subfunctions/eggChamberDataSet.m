
classdef eggChamberDataSet < handle

    properties (GetAccess = 'public', SetAccess = 'public')
        
        % table holding the actual data
        t = table();

        % exhaustive dataset before removing useless variables and metrics
        % left in case in-depth QC is needed.
        fullT = table();
        
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

        % ******************** data I/O properties

        %  name of the subfolder where the egg chamber data is stored
        eggChamberCsvFolderName = 'eggChamberCSV';
    
        % name of the subfolder where the 2D segmentations of the egg chambers are stored
        % should hold a file called eggChamberStages.csv - it is optional.
        eggChamberSegFolderName = 'eggChamberSEG';
        eggChamberSegFileName = 'eggChamberStages.csv';
    
        % file names (of name base in the case of channel dependent file names)
        eggChamberWiseFileNames = {'allNucGeom.csv';...
            '_allNucInt.csv';...
            '_wholeImgInt.csv';...
            '_eggChamberInt.csv'};
    
        % whether each file name of the previous list is channel dependent or not
        eggChamberWiseChannelDependent = logical([0;...
            1;...
            1;...
            1;]);
    
        % whether each file name of the previous list contains a single row (whole image data) or
        % one row per nucleus
        eggChamberWiseSingleRow = logical([0;...
            0;...
            1;...
            1;]);
        
        % prefix to add in front of each variable name when loading the
        % corresponding table
        eggChamberWisePrefix ={'nuc';...
            'nuc';...
            'wholeImg';...
            'eggChamber'};
    
        % prefix to add at the end of each variable name when loading the
        % corresponding table
        eggChamberWiseSuffix ={'';...
            '_raw';...
            '_raw';...
            '_raw'};

        % ******************** Nuclei Table variables properties

        % the variables (i.e. columns) which store the sample and condition ID
        sampleVarList = {'condIdx', 'sampleIdx','nucInputFileName'};
        
        % prefixes 
        % in a variable name, the prefix indicates what structure/ROI type the metric is applied to
        nucPrefix = 'nuc';
        wholeImgPrefix = 'wholeImg';
        eggChamberPrefix = 'eggChamber';
        clusterPrefix = 'clust';
        prefixList;
    
        % list of the variable basenames that pertain to geometry - note
        % that the eggChamberID and eggChamberStage fall into that category
        geomVarsBaseNameList = {'Label', 'Volume', 'SurfaceArea',...
            'MeanBreadth', 'Sphericity', 'EulerNumber', ...
        'Box_X_Min', 'Box_X_Max', 'Box_Y_Min', 'Box_Y_Max', ...
        'Box_Z_Min', 'Box_Z_Max', 'Centroid_X', 'Centroid_Y', 'Centroid_Z', ...
        'Dist',...
        'eggChamberID','eggChamberStage'}; 
        % make sure that no entry in geomVarsBaseNameList is also part of the name of one of the entries in channelVarsBaseNameList
    
        % list of the variable basenames that pertain to an intensity metric
        channelVarsBaseNameList = {'Mean', 'StdDev', 'Max', 'Min','Median', 'Mode', 'Skewness', 'Kurtosis','Volume','NumberOfVoxels'}; 
        % the fact that 'Mean' is also part of a variable ('MeanBreadth') in geomVarsBaseNameList shouldnt be an issue
        % note that Volume is both an intensity and a geometry metric,
        % annoyingly

        % name of variables to remove in streamlined table
        geomVarsToRemoveInStreamLinedNucleiTable = {'SurfaceArea','MeanBreadth', 'Sphericity', 'EulerNumber', ...
        'Box_X_Min', 'Box_X_Max', 'Box_Y_Min', 'Box_Y_Max', ...
        'Box_Z_Min', 'Box_Z_Max', 'Centroid_X', 'Centroid_Y', 'Centroid_Z', ...
        'Dist'};
        channelVarsToRemoveInStreamLinedNucleiTable = {'Mode', 'Skewness', 'Kurtosis','StdDev','Volume','NumberOfVoxels'};
    
        % suffixes
        % in an intensity variable name, the suffix indicates the processing applied to
        % the intensity
        rawSuffix = 'raw';
        plasmCorrSuffix = 'plasmCorr';
        eggChamberCorrSuffix = 'eggChamberCorr';
        wholeImgCorrSuffix = 'wholeImgCorr';
        suffixList;
    
        % variable names are built as follows:
            % sample ID variables: the full name is in the list sampleVarList
    
            % geometry variables: <prefix><geomVarBaseName>
    
            % channel variables:
            % <prefix>C<channelIdx>_<channelVarsBaseNameList>_<suffix>
        
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

            obj.prefixList = {obj.nucPrefix,obj.wholeImgPrefix,obj.eggChamberPrefix,obj.clusterPrefix};
            obj.suffixList = {obj.rawSuffix,obj.plasmCorrSuffix,obj.eggChamberCorrSuffix,obj.wholeImgCorrSuffix};
            obj.backgroundIntensityPrefixList = {obj.wholeImgPrefix,obj.eggChamberPrefix};

        end
        
        %% load data from all conditions into data table
        function loadAllEggChamberNucleiData(obj)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames'); % mute unneeded warnings
            removeNeighbors = 1; % exclude from the data import the variables with neighbor in their name - we do not use them.
            obj.t = [];
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    disp(['Loading data from condition ',obj.conditionNames{i},...
                        ', sample ',obj.sampleNames{i}{j},' ...']);
                    obj.t = obj.combineEcTables(obj.t, loadEggChamberData(obj,i,j,removeNeighbors));
                end
            end
            obj.fullT = obj.t; % fullT backup copy of the exhaustive imported data table
            obj.getChannelList;
            obj.getEggChamberIDs;
            disp('done.');
        end
        
        

        %% remove unlikely to be used variables from Nuclei table
        function streamLineNucleiTable(obj)
            disp('Streamlining Nuclei table...');
            [c,nChannels] = obj.getChannelList;
            
            if ismember('eggChamberID',obj.t.Properties.VariableNames) ...
                    && ismember('eggChamberStage',obj.t.Properties.VariableNames) ...
                    && obj.eggChamberSegChannel ~= 0
                removeEggChamberSegChannelVars = 1;
            else
                removeEggChamberSegChannelVars = 0;
            end
            
            % compile list of variables to remove
            varsToRemoveList = {};
            for i=1:numel(obj.prefixList)
                for n=1:numel(obj.geomVarsToRemoveInStreamLinedNucleiTable)
                    newVar = buildVarName(obj,...
                        obj.prefixList{i},0,obj.geomVarsToRemoveInStreamLinedNucleiTable{n},'','geom');
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
                            curVarList = obj.channelVarsToRemoveInStreamLinedNucleiTable;
                        end

                        for n=1:numel(curVarList)
                            newVar = buildVarName(obj,...
                                obj.prefixList{i},c(j),...
                                curVarList{n},obj.suffixList{k},'channel');
                            varsToRemoveList = [varsToRemoveList; newVar];
                        end
                    end
                end
            end

            % remove variable names that arent present
            varsToRemoveList = varsToRemoveList(...
                ismember(varsToRemoveList,obj.t.Properties.VariableNames));

            obj.t = removevars(obj.t,varsToRemoveList);

            % reorder variables
            obj.reOrderNucleiVariables();

            % sort rows
            obj.sortNucleiRowsByEggChamber();

            disp('Done.');
        end

        %% perform background subtraction on nuclei intensity values
         function backgroundCorrectNucIntensity(obj)
            
            obj.t = bgCorrTable(obj.t);
            obj.fullT = bgCorrTable(obj.fullT);

            function tOut = bgCorrTable(tIn)
                tOut = tIn;

                [c,nChannels] = obj.getChannelList;
                
                % collect the background columns to subtract
                varToSubtract = cell(...
                    nChannels,numel(obj.backgroundIntensityPrefixList));
    
                for i=1:nChannels
                    for j=1:numel(obj.backgroundIntensityPrefixList)
                        
                        varToSubtract{i,j} = [obj.backgroundIntensityPrefixList{j},...
                            'C',num2str(c(i)),'_',obj.metricToSubtract,'_',obj.rawSuffix];
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
                        repmat({['C',num2str(c(i)),'_']}, size(varList)),...
                        'UniformOutput',0) );
                
                    curVarList = varList(idxC);
                    for j=1:numel(obj.backgroundIntensityPrefixList)
                        for k =1:numel(curVarList)
                            newVarName = strrep(curVarList{k},...
                                ['_',obj.rawSuffix],...
                                ['_',obj.backgroundIntensityPrefixList{j},'Corr']);
    
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
            eggChamberID = [];
            eggChamberStage = [];
            numberOfNucleiInEggChamber = [];
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    for k=1:obj.eggChamberNumber{i}(j)
                        conditionName = ...
                            [conditionName; obj.conditionNames{i}];
                        sampleName = [sampleName; obj.sampleNames{i}{j}];
                        eggChamberID = ...
                            [eggChamberID; obj.eggChamberIDs{i}{j}(k)];
                        eggChamberStage = ...
                            [eggChamberStage; ...
                            obj.eggChamberStages{i}{j}(k)];
                        numberOfNucleiInEggChamber = ...
                            [numberOfNucleiInEggChamber; ...
                            obj.eggChamberNumNucPerEC{i}{j}(k)];
                    end
                end
            end
            sumT = table(conditionName,sampleName,eggChamberID,...
                eggChamberStage,numberOfNucleiInEggChamber);
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
                    x = obj.t.(varName)(...
                        obj.t.condIdx ==obj.condIndices(j) ...
                        & obj.t.sampleIdx == s(k));
            
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
                        curYPlot = obj.t.(varName)(...
                            obj.t.condIdx ==obj.condIndices(j) ...
                            & obj.t.sampleIdx == s(k))';
            
                    elseif nNuclei == 1
                        curXPlot = xSampleVals{j}(k);
                        curYPlot = obj.t.(varName)(...
                            obj.t.condIdx ==obj.condIndices(j) ...
                            & obj.t.sampleIdx == s(k))';
            
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
            if ~ismember( varName, obj.t.Properties.VariableNames)
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
                            x = obj.t.(varName)(...
                                obj.t.condIdx ==obj.condIndices(i) ...
                                & obj.t.sampleIdx == s(j) ...
                                & obj.t.eggChamberID == obj.eggChamberIDs{i}{j}(k));
                    
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
                                curYPlot = obj.t.(varName)(...
                                    obj.t.condIdx ==obj.condIndices(i) ...
                                    & obj.t.sampleIdx == s(j)...
                                    & obj.t.eggChamberID == obj.eggChamberIDs{i}{j}(k))';
                    
                            elseif nNuclei == 1
                                curXPlot = xEggChamberVals{i}{j}(k);
                                curYPlot = obj.t.(varName)(...
                                    obj.t.condIdx ==obj.condIndices(i) ...
                                    & obj.t.sampleIdx == s(j)...
                                    & obj.t.eggChamberID == obj.eggChamberIDs{i}{j}(k))';
                    
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
   
    end

    methods(Access = 'private')


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
        function [sampleIdx,conditionIdx] = getSampleIdx(obj,cName,sName)
            conditionIdx = find(ismember( obj.conditionNames,cName));
            sampleIdx = find(ismember( obj.sampleNames{conditionIdx},sName));
        end
    
        %% get the list of expected files to be generated by a complete
        % analysis of the eggChamber data (only channel-independent files).
        function fNames = expectedEggChamberAnalysisFilesChannelIndependentFiles(obj)
            fNames = obj.eggChamberWiseFileNames(...
                ~obj.eggChamberWiseChannelDependent);
        end
    
        %% get the list of expected files to be generated by the
        % analysis of the eggChamber data for a specific channel .
        function fNames = expectedEggChamberAnalysisFilesChannelFiles(obj,channelIdx)
           
            nameList = obj.eggChamberWiseFileNames(...
                obj.eggChamberWiseChannelDependent);
    
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
        function analysisDir = getEggChamberAnalysisDir(obj,conditionIdx,sampleIdx)
            if isempty(obj.conditionNames)
                disp('condition folder appears empty.');
                analysisDir = [];
                return ;
            end
            analysisDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sampleIdx},...
                obj.eggChamberCsvFolderName);
        end
    
        function eggSegDir = getEggChamberSegDir(obj,conditionIdx,sampleIdx)
            if isempty(obj.conditionNames)
                disp('condition appears empty.');
                eggSegDir = [];
                return ;
            end
            eggSegDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sampleIdx},...
                obj.eggChamberSegFolderName);
        end
        
        %% check whether all files expected from the analysis of the 
        % eggchamber data are present in the analysis folder
        function [channelsFound,eggChamberSegFound] = areEggChamberAnalysisFilesPresent(obj,...
                                conditionIdx,sampleIdx)
    
            % arbitrary max number of channels to check
            maxChannels = 10;
    
            % check that the analysis dir is present
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sampleIdx);
            
            if ~exist(analysisDir, 'dir')
                channelsFound = 0;
                eggChamberSegFound = 0;
                return;
            end
    
            % check that the egg chamber seg dir is present
            eggSegDir = obj.getEggChamberSegDir(conditionIdx,sampleIdx);
            
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
            for i=1:maxChannels
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
        % conditionIDx, sampleIdx relative to conditionNames and
        % sampleNames cell arrays. 
        % removeNeighbors is a flag - set to 1 to remove variables that contain the word neighbor
        % (these are generated by Fiji but tend to not be useful for our
        % analyses).
        function t = loadEggChamberData(obj,conditionIdx,sampleIdx,removeNeighbors)
            
            % check that data is present
            [c,eggSegFound] = obj.areEggChamberAnalysisFilesPresent(conditionIdx,sampleIdx);
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
    
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sampleIdx);
            
            tS = []; % full table place holder for data that is single row 
            % - e.g. whole image intensity, whole egg chamber intensity.
            t = []; % full table place holder for data that is mutli row
            % i.e. nucleus-specific 
            % loop through files, load each as a table curT and append the 
            % variables curT holds as new columns into t (or tS if single
            % row)
            for i=1:numel(obj.eggChamberWiseFileNames)
                if ~obj.eggChamberWiseChannelDependent(i)
                    % load table
                    curT = readtable( fullfile(analysisDir,...
                        obj.eggChamberWiseFileNames{i}) );
                    
                    % remove useless variable name (if present)
                    if ismember('Var1',curT.Properties.VariableNames)
                        curT = removevars( curT,'Var1');
                    end
    
                    % add prefix and suffix to variable name
                    for k=1:numel(curT.Properties.VariableNames)
                        if ~strcmp(curT.Properties.VariableNames{k},'Label')
                            curT.Properties.VariableNames{k} = ...
                                [obj.eggChamberWisePrefix{i},...
                                curT.Properties.VariableNames{k},...
                                obj.eggChamberWiseSuffix{i}];
                        end
                    end
                    
                    % append to full table
                    if obj.eggChamberWiseSingleRow(i)
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
                            obj.eggChamberWiseFileNames{i}]) );
                        
                        % remove useless variable name (if present)
                        if ismember('Var1',curT.Properties.VariableNames)
                            curT = removevars( curT,'Var1');
                        end
    
                        % add prefix and suffix to variable name
                        for k=1:numel(curT.Properties.VariableNames)
                            if ~strcmp(curT.Properties.VariableNames{k},'Label')
                                curT.Properties.VariableNames{k} = ...
                                    [obj.eggChamberWisePrefix{i},...
                                    'C',num2str(j),'_',...
                                    curT.Properties.VariableNames{k},...
                                    obj.eggChamberWiseSuffix{i}];
                            end
                        end
    
                        % append to full table
                        if obj.eggChamberWiseSingleRow(i)
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
    
            % rename Label variable to nucLabel
            t =renamevars(t,{'Label'},{'nucLabel'});
    
            % add condition and sample index variables
            cIdx = repmat(conditionIdx,size(t,1),1);
            sIdx = repmat(sampleIdx,size(t,1),1);
            t = addvars(t,cIdx,sIdx,...
                'NewVariableNames',{'condIdx','sampleIdx'},...
                'Before',1);
    
            % if egg chamber stage file is present, add a variable for the egg
            % chamber stage.
            if eggSegFound && (obj.eggChamberSegChannel ~= 0)
                % get the name of the variable holding the median nucleus
                % intensity in the channel that holds the egg chamber
                % segmenentation ID (the median int will be used as the ID).
                eggChamberIDVariable = [obj.eggChamberWisePrefix{2},'C',...
                    num2str(obj.eggChamberSegChannel),'_Median',obj.eggChamberWiseSuffix{2}];
                
                % convert the egg chamber ID into and egg chamber stage based
                % on the key saved in the file.
                ecStage = obj.getEggChamberStageData(t.(eggChamberIDVariable),conditionIdx,sampleIdx);
                
                % copy the eggChamber variable with ean easy to interpret file
                % name
                t = addvars(t,t.(eggChamberIDVariable),...
                'NewVariableNames',{'eggChamberID'});
    
                % add egg Chamber stage variable to table
                t = addvars(t,ecStage,...
                'NewVariableNames',{'eggChamberStage'});
            end
    
            % remove variables that contain the word neighbor
            if removeNeighbors
                idx = cellfun(@contains,t.Properties.VariableNames,...
                    repmat( {'Neighbors'},size(t.Properties.VariableNames) ),'UniformOutput',0);
                idx = cell2mat(idx);
                t = t(:,~idx);
            end
        end
        
        %% concatenate tables vertically. If variables are missing in one of the
        % tables, it is added with intialized defaulted values to zero or ''.
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
                        newVar = zeros(size(t2,1),1);
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
                        newVar = zeros(size(t1,1),1);
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
        function ecStage = getEggChamberStageData(obj, eggChamberID,conditionIdx,sampleIdx)
            % load cell stage key file
            eggSegDir = obj.getEggChamberSegDir(conditionIdx,sampleIdx);
            ecStageFile = fullfile(eggSegDir,obj.eggChamberSegFileName);
    
            % reads the egg chamber stage key file as a 2-column table w
            % variable names eggChamberID and eggChamberStage.
            tEC = readtable(ecStageFile);
    
            % convert key file to Containers.map object.
            ecMap = containers.Map(tEC.eggChamberID,tEC.eggChamberStage);
            ecStage = zeros(numel(eggChamberID),1);
            ecVals = unique(eggChamberID);
            for i=1:numel(ecVals)
                if ~isKey(ecMap,ecVals(i))
                    ecStage(eggChamberID==ecVals(i)) = 0;
                else
                    ecStage(eggChamberID==ecVals(i)) = cell2mat(values(ecMap,{ecVals(i)}));
                end
            end
        end



        %%
        function reOrderNucleiVariables(obj)
            % order should go:
            % input file name
            % condIdx
            % sampleIdx
            % eggChamberID (optional)
            % eggchamber stage (optional)
            % nucLabel
            % nucVolume 
            
            vStart = obj.t.Properties.VariableNames{1};
            obj.t = movevars(obj.t,'nucInputFileName','Before',vStart);
            obj.t = movevars(obj.t,'condIdx','After','nucInputFileName');
            obj.t = movevars(obj.t,'sampleIdx','After','condIdx');
            if ismember('eggChamberID',obj.t.Properties.VariableNames)
                obj.t = movevars(obj.t,'eggChamberID','After','sampleIdx');
            end

            if ismember('eggChamberStage',obj.t.Properties.VariableNames)
                obj.t = movevars(obj.t,'eggChamberStage','After','eggChamberID');
                obj.t = movevars(obj.t,'nucLabel','After','eggChamberStage');
            else
                if ismember('eggChamberID',obj.t.Properties.VariableNames)
                    obj.t = movevars(obj.t,'nucLabel','After','eggChamberID');
                else
                    obj.t = movevars(obj.t,'nucLabel','After','sampleIdx');
                end
            end
            obj.t = movevars(obj.t,'nucVolume','After','nucLabel');
        end

        %%
        function sortNucleiRowsByEggChamber(obj)
            if ismember('eggChamberID',obj.t.Properties.VariableNames)
                obj.t = sortrows(obj.t,{'condIdx','sampleIdx','eggChamberID','nucLabel'});
            else
                disp('Cannot sort table by egg chamber, eggChamberID column missing. Skipping.');
            end

        end

        %% find the number of channels in the table
        function [channelIndices,nChannels] = getChannelList(obj)
            allVars = obj.t.Properties.VariableNames;
            channelIndices = [];
            for i=1:numel(allVars)
                v = obj.getVarType(allVars{i});
                if strcmp(v,'channel')
                    [~,c,~,~] = obj.splitVarName(allVars{i});
                    channelIndices = [channelIndices,c];
                end
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
                    p = varName(1:k-1); % p should be something like nucC1_
                    for i=1:numel(obj.prefixList)
                        if contains(p,obj.prefixList{i})
                            prefix = obj.prefixList{i};
                            % now that we found the prefix in our string (say 'nuc' in 'nucC1_'), we
                            % collect the channel ID by starting two
                            % characters after the end of the prefix, and
                            % removing the last character ('_'):
                            channel = str2double(p(length(prefix)+2:end-1));
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

            switch varType
                case 'sample'
                    varName = baseName;
                case 'geom'
                    if ~ismember(prefix,obj.prefixList)
                        print('could not build variable name, prefix ',prefix,' not recognized');
                        return
                    else
                        varName = [prefix,baseName];
                    end
                case 'channel'
                    if ~ismember(prefix,obj.prefixList)
                        print('could not build variable name, prefix ',prefix,' not recognized');
                        return
                    else
                        obj.channelList = obj.getChannelList;
                        if ~ismember(channel,obj.channelList)
                            print('could not build variable name, channel ',num2str(channel),' not recognized');
                            return
                        end
    
                        if ~ismember(suffix,obj.suffixList)
                            print('could not build variable name, suffix ',suffix,' not recognized');
                            return
                        end
    
                        varName = [prefix,'C',num2str(channel),'_',baseName,'_',suffix];
                    end
            end        
        end
        
         %% find whether a given variable (varName) is a sample variable, a
        % geometry variable or a channel variable. ALso outputs the
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
        % ecN, ecID and ecStage are defaulted to {} if the eggChamberID
            % variable is absent from the data table.
        function  [ecID,ecN,ecStage,ecNumNucPerEC] = getEggChamberIDs(obj)
            % check whether the eggChamberID variable is present. If not, return
            % empty cell arrays.
            if ~ismember('eggChamberID',obj.t.Properties.VariableNames)
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
                        unique(obj.t.eggChamberID(...
                        obj.t.condIdx == i & obj.t.sampleIdx == j));
                    ecN{i}(j) = numel(ecID{i}{j});

                    curT = obj.t( obj.t.condIdx == i & obj.t.sampleIdx == j,:);
                    for k=1:numel(ecID{i}{j})
                        ecNumNucPerEC{i}{j}(k) = sum(curT.eggChamberID == ecID{i}{j}(k));
                        if ismember('eggChamberStage',obj.t.Properties.VariableNames)
                            ecStage{i}{j}(k) = ...
                            unique(...
                            curT.eggChamberStage(curT.eggChamberID == ecID{i}{j}(k)));
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
    end
end
