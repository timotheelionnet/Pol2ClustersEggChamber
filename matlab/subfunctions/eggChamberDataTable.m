
classdef eggChamberDataTable < handle

    properties (GetAccess = 'public', SetAccess = 'public')
        
        % the actual data table
        t;
        
        % list of channels
        channelList;

        % list of condition indices
        condIndices;

        % list of sample indices and the number of samples per condition
        sampleIndices;
        nSamples;
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
    
        % list of the variable basenames that pertain to geometry
        geomVarsBaseNameList = {'Label', 'Volume', 'SurfaceArea',...
            'MeanBreadth', 'Sphericity', 'EulerNumber', ...
        'Box_X_Min', 'Box_X_Max', 'Box_Y_Min', 'Box_Y_Max', ...
        'Box_Z_Min', 'Box_Z_Max', 'Centroid_X', 'Centroid_Y', 'Centroid_Z', ...
        'Dist'}; 
        % make sure that no entry in geomVarsBaseNameList is also part of the name of one of the entries in channelVarsBaseNameList
    
        % list of the variable basenames that pertain to an intensity metric
        channelVarsBaseNameList = {'Mean', 'StdDev', 'Max', 'Min','Median', 'Mode', 'Skewness', 'Kurtosis','StdDev','NumberOfVoxels'}; 
        % the fact that 'Mean' is also part of a variable ('MeanBreadth') in geomVarsBaseNameList shouldnt be an issue
    
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
        % used to set the spacing between cloud of points    
        spacingUnit = 1;

        % how much separation there is between conditions as a function of the separation between samples.
        condSeparator = 1.5; 
        
        % free space between samples
        freeSpaceBetweenSamples = 0.2;
        
    end

    methods
        %% initialize object
        function obj = eggChamberDataTable(t)
            obj.t = t;
            obj.prefixList = {obj.nucPrefix,obj.wholeImgPrefix,obj.eggChamberPrefix,obj.clusterPrefix};
            obj.suffixList = {obj.rawSuffix,obj.plasmCorrSuffix,obj.eggChamberCorrSuffix,obj.wholeImgCorrSuffix};
            obj.backgroundIntensityPrefixList = {obj.wholeImgPrefix,obj.eggChamberPrefix};

            obj.getChannelList;
            obj.getConditions;
            obj.getSamples;
        end

        %% perform background subtraction on nuclei intensity values
         function obj = backgroundCorrectNucIntensity(obj)
            t2 = obj.t;
               
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
            varList = obj.t.Properties.VariableNames;
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
            
            varList = obj.t.Properties.VariableNames(idx);
            
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

                        newVar = obj.t.(curVarList{k}) - obj.t.(varToSubtract{i,j});
                        t2 = addvars(t2,newVar,'NewVariableNames',newVarName);
                    end
                end
            end
            obj.t = t2;
        end

        %% scatter plot a metric by sample
        function scatterPlotMetricBySample(obj,prefix,channel,baseName,suffix)

            varName = obj.buildVarName(prefix,channel,baseName,suffix);
            
            % make sure conditions/samples numbers are up to date
            obj.getConditions;
            obj.getSamples;

            % build figure
            figure('Name',strrep(varName,'_',' '));
            hold;
            
            % collect the x values to plot each condition/sample at
            [xSampleVals,xSampleIDs] = obj.getSampleXValues(...
                obj.condIndices,obj.nSamples,obj.spacingUnit);

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
                    x = obj.t.(varName)(obj.t.condIdx ==obj.condIndices(j) & obj.t.sampleIdx == s(k));
            
                    % generate slightly offset x coordinates for each nucleus,
                    % centered around the sample X
                    nNuclei = size(x,1);
                    if nNuclei >1
                        % spacing between nuclei
                        nucSpacing = ...
                            obj.spacingUnit*(1-2*obj.freeSpaceBetweenSamples)/(nNuclei-1);
            
                        % x coordinate for each nucleus of current condition/sample
                        curXPlot = xSampleVals(j,k) - floor(nNuclei/2)*nucSpacing ...
                            + (0:(nNuclei-1))*nucSpacing;
            
                        % y coordinate for each nucleus of current condition/sample
                        curYPlot = obj.t.(varName)(obj.t.condIdx ==obj.condIndices(j) & obj.t.sampleIdx == s(k))';
            
                    elseif nNuclei == 1
                        curXPlot = xSampleVals(j,k);
                        curYPlot = obj.t.(varName)(obj.t.condIdx ==obj.condIndices(j) & obj.t.sampleIdx == s(k))';
            
                    elseif nNuclei == 0
                        curXPlot = [];
                        curYPlot = [];
                    end
            
                    % append coordinates of current condition/sample to global list
                    xPlot = [xPlot,curXPlot];
                    yPlot = [yPlot,curYPlot];
            
                    xSampleValsVec = [xSampleValsVec,xSampleVals(j,k)];
                    xSampleIDsVec = [xSampleIDsVec,xSampleIDs{j,k}];
                end
            end
            plot(xPlot,yPlot,'o');
            xticks(xSampleValsVec);
            xticklabels(xSampleIDsVec);
            xtickangle(45);
            ylabel(baseName);

        end

        %% get the list of indices for the conditions in the table (conditionIndices), and
        % their number (nConditions)
        function  [condIdxList,nConds] = getConditions(obj)
            obj.condIndices = unique(obj.t.condIdx);
            nConds = numel(obj.condIndices);
            condIdxList = obj.condIndices;
        end
        
        %% get the list of sample indices for each condition in the table 
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
        function varName = buildVarName(obj,prefix,channel,baseName,suffix)

            varName = [];
            
            if ismember(baseName,obj.sampleVarList) 
                varName = baseName;
                return
            end

            if ismember(baseName,obj.geomVarsBaseNameList)
                if ~ismember(prefix,obj.prefixList)
                    print('could not build variable name, prefix ',prefix,' not recognized');
                    return
                else
                    varName = [prefix,baseName];
                end
            end

            if ismember(baseName,obj.channelVarsBaseNameList)
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
        
        %% generate x coordinates by sample/condition for scatter plots
        function [xSampleVals, xSampleIDs] = getSampleXValues(obj,conditionNames,nSamples,varargin)
        % generates a series of X coordinates for each sample, separated by
        % condition xSampleVals(i,j) is the x coordinate for condition i, sample j 
        % and matching tick names xSampleIDs{i,j}
        
        % conditionNames: 1 x n cell array containing the names of the various
            % conditions.
        % nSamples 1 x n array containing the number of samples per each condition.
        % optional argument spacingUnit is the unit spacing between two samples - its value is irrelevant 
            % unless you want to use the X values in a plot where other things are present and need 
            % to coordinate the X values. default if not populated: 1.
        
        
        if numel(varargin) ~= 0
            obj.spacingUnit = varargin{1};
        end

        if isa(conditionNames,'numeric')
            conditionNames = num2cell(conditionNames);
            conditionNames = cellfun(@num2str,conditionNames,'UniformOutput',0);
        end
        
        nConditions = numel(conditionNames);
        
        xSampleVals = zeros(nConditions,max(nSamples));
        xSampleIDs = cell(nConditions,max(nSamples));
        curX = 0;
        for i=1:nConditions
            if i>1
                curX = curX + obj.condSeparator*obj.spacingUnit;
            end
            
            xSampleVals(i,1:nSamples(i)) = curX + (1:nSamples(i))*obj.spacingUnit;
            curX = curX + nSamples(i)*obj.spacingUnit;
        
            for j=1:nSamples(i)
                xSampleIDs{i,j} = [conditionNames{i},' sample',num2str(j)];
            end
        end


end
    end
end
