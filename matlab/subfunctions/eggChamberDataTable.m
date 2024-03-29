
classdef eggChamberDataTable < handle

    properties (GetAccess = 'public', SetAccess = 'public')
        
        % table holding the actual data
        t;

        % exhaustive dataset before removing useless variables and metrics
        % left in case in-depth QC is needed.
        fullT;
        
        % list of channels, e.g. [1,2,3,4]
        channelList;

        % list of condition indices (num array) and corresponding folder names (cell array) 
        % condIndices is a nConditions x 1 numerical array which holds the
            % index of each condition  (usually will be 1:nConditions in a vertical vector)
            % e.g. [1;2]
        condIndices;
        % conditionNames is a nConditions x 1 cell array where the ith
            % entry is the name of the ith condition
            % e.g. {'Ctrl';
            %       'TRI';}
        conditionNames;

        % list of sample indices (num array) and the number of samples per condition
        % nSamples is an nConditions x 1 numerical array where ith entry 
            % is the number of samples in condition i).
            % e.g. [3;3]
        nSamples;
        % sampleIndices is an nConditions x 1 cell array where ith entry is
            % an nSamples(i) x 1 numeric array holding the index of each
            % sample (usually will be 1:nSamples(i) in a vertical vector)
            % e.g. sampleIndices{2,1} = [1;2;3]
        sampleIndices;
        % sampleNames is an nConditions x 1 cell array where ith entry is
            % a cell array (nSamples(i) x1) where the jth entry is the name
            % of the sample
            % e.g. sampleNames{2,1} = {'9-TRI-646,MPM2-488,Ser5ph-Cy3-1001zCorr';
            %                          '9-TRI-646,MPM2-488,Ser5ph-Cy3-1002zCorr';
            %                          '9-TRI-646,MPM2-488,Ser5ph-Cy3-1zCorr';}
        sampleNames; 

        % input Folder path, e.g. '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img'
        inFolder;

        % index of the data channel that holds the egg chamber ID (scalar; set to zero if no
            % such info available).
            % e.g. 5
        eggChamberSegChannel;
        
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
        
        %% nucleus intensity background subtraction settings
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

        %% plotting settings    
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

        % name of variables to remove in streamlined table
        geomVarsToRemoveInStreamLinedNucleiTable = {'SurfaceArea','MeanBreadth', 'Sphericity', 'EulerNumber', ...
        'Box_X_Min', 'Box_X_Max', 'Box_Y_Min', 'Box_Y_Max', ...
        'Box_Z_Min', 'Box_Z_Max', 'Centroid_X', 'Centroid_Y', 'Centroid_Z', ...
        'Dist'};
        channelVarsToRemoveInStreamLinedNucleiTable = {'Mode', 'Skewness', 'Kurtosis','StdDev','Volume','NumberOfVoxels'};
    end

    methods (Access = 'public')
        %% initialize object
        % t is the data table, ec is an object of eggChamberDataFolder
        % class, which containts the location
        function obj = eggChamberDataTable(t,ec)
            obj.t = t;
            obj.fullT = t; % fullT backup copy of the exhaustive imported data table
            obj.prefixList = {obj.nucPrefix,obj.wholeImgPrefix,obj.eggChamberPrefix,obj.clusterPrefix};
            obj.suffixList = {obj.rawSuffix,obj.plasmCorrSuffix,obj.eggChamberCorrSuffix,obj.wholeImgCorrSuffix};
            obj.backgroundIntensityPrefixList = {obj.wholeImgPrefix,obj.eggChamberPrefix};

            obj.getChannelList;
            obj.getConditions;
            obj.getSamples;
            
            obj.inFolder = ec.inFolder;
            obj.conditionNames = ec.conditionNames;
            obj.sampleNames = ec.sampleNames;
            obj.eggChamberSegChannel = ec.eggChamberSegChannel;
            obj.getEggChamberIDs;
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
         function obj = backgroundCorrectNucIntensity(obj)
            
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
            [~,nConds] = obj.getConditions;
            for i=1:nConds
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
            
            % make sure conditions/samples numbers are up to date
            obj.getConditions;
            obj.getSamples;

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

            % make sure conditions/samples numbers are up to date
            obj.getConditions;
            obj.getSamples;

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

        %% get the list of indices for the conditions in the table (conditionIndices), and
        % their number (nConditions)
        function  [condIdxList,nConds] = getConditions(obj)
            obj.condIndices = unique(obj.t.condIdx);
            nConds = numel(obj.condIndices);
            condIdxList = obj.condIndices;
        end      
       
        %% generate x coordinates by sample/condition for scatter plots
        function [xSampleVals, xSampleIDs] = getSampleXValuesBySample(obj)
            % subfunction used to generate x coordinates for plotting,
            % grouped by sample, then condition.
            % xSampleVals{i}(j) is the x coordinate for condition i, sample j 
            % and matching tick names xSampleIDs{i}{j}
            
            nConditions = numel(obj.conditionNames);
            
            xSampleVals = cell(nConditions,1);
            xSampleIDs = cell(nConditions,1);
            curX = 0;
            for i=1:nConditions
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
            nConditions = numel(obj.conditionNames);
            xEggChamberVals = cell(nConditions,1);
            xEggChamberIDs = cell(nConditions,1);
            xEggChamberStages = cell(nConditions,1);
            curX = 0;
            for i=1:nConditions
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
        
        %% get the list of sample indices for each condition in the table 
        % outputs sI which is the updated obj.sampleIndices 
            % (nConditions x 1 cell array where ith entry is a numeric array 1:nSamples(i))
            % outputs nS which is the updated obj.nSamples (nConditions x 1
            % numerical array where ith entry is the number of samples in
            % condition i).
        function  [sI,nS] = getSamples(obj)
            % make sure list of conditions is up to date
            obj.condIndices = obj.getConditions;
            
            obj.sampleIndices = cell(numel(obj.condIndices),1);
            obj.nSamples = zeros(numel(obj.condIndices),1);
            for i=1:numel(obj.condIndices)
                obj.sampleIndices{i,1} = unique(obj.t.sampleIdx(obj.t.condIdx == obj.condIndices(i))); 
                obj.nSamples(i,1) = numel(obj.sampleIndices{i,1});
            end
            sI = obj.sampleIndices;
            nS = obj.nSamples;
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

            % make sure list of conditions and samples is up to date
            obj.condIndices = obj.getConditions;
            obj.getSamples;

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
