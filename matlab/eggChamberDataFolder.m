classdef eggChamberDataFolder < handle
% from the below expected folder architecture, collects conditions and
% samples

% -input folder (arbitrary name)
%   - condition1 (arbitrary name)
%       - eggChamber1 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%       - eggChamber2 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%   - condition2
%       - eggChamber1 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%       - eggChamber2 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 

properties (GetAccess = 'public', SetAccess = 'public')
    inFolder;
    conditionNames;
    nConditions;
    sampleNames;
    nSamples;
end

properties (GetAccess = 'private', SetAccess = 'private')

    %  name of the subfolder where the egg chamber data is stored
    eggChamberCsvFolderName = 'eggChamberCSV';

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
end

methods
        %% initialize object
        function obj = eggChamberDataFolder(inputFolder)
            obj.inFolder = inputFolder;
            obj.collectConditionsAndSamples;
        end

        %% collect samples and conditions
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

            % for each condition, collect samples, i.e. list of subfolders
            % within the condition folder:
            for i=1:obj.nConditions
                d = dir(fullfile(obj.inFolder, obj.conditionNames{i}));
                dFolders = d([d(:).isdir]);
                dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
                obj.sampleNames{i,1} = {dFolders(:).name};
                obj.sampleNames{i,1} = reshape(obj.sampleNames{i},...
                    numel(obj.sampleNames{i}),1);
                obj.nSamples(i,1) = numel(dFolders);
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
                disp('eggChmaber appears empty.');
                analysisDir = [];
                return ;
            end
            analysisDir = fullfile(obj.inFolder,...
                obj.conditionNames{conditionIdx},...
                obj.sampleNames{conditionIdx}{sampleIdx},...
                obj.eggChamberCsvFolderName);
        end
        
        %% check whether all files expected from the analysis of the 
        % eggchamber data are present in the analysis folder
        function channelsFound = areEggChamberAnalysisFilesPresent(obj,...
                                conditionIdx,sampleIdx)

            % arbitrary max number of channels to check
            maxChannels = 10;

            % check that the analysis dir is present
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sampleIdx);
            
            if ~exist(analysisDir, 'dir')
                channelsFound = 0;
                return;
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
        
        %%
        function t = loadEggChamberData(obj,conditionIdx,sampleIdx,removeNeighbors)
            
            % check that data is present
            c = obj.areEggChamberAnalysisFilesPresent(conditionIdx,sampleIdx);
            if c == 0
                nr = [];
                return
            end
            
            analysisDir = obj.getEggChamberAnalysisDir(conditionIdx,sampleIdx);
            
            tS = [];
            t = [];
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
                    if obj. eggChamberWiseSingleRow(i)
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
                        if obj. eggChamberWiseSingleRow(i)
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

            % remove variables that contain the word neighbor
            if removeNeighbors
                idx = cellfun(@contains,t.Properties.VariableNames,...
                    repmat( {'Neighbors'},size(t.Properties.VariableNames) ),'UniformOutput',0);
                idx = cell2mat(idx);
                t = t(:,~idx);
            end
        end

        %%
        function t = loadAllEggChamberData(obj)
            warning('off','MATLAB:table:ModifiedAndSavedVarnames');
            removeNeighbors = 1;
            t = [];
            for i=1:obj.nConditions
                for j=1:obj.nSamples(i)
                    t = [t; loadEggChamberData(obj,i,j,removeNeighbors)];
                end
            end
        end

end

end