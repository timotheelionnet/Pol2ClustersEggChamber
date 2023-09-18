% enter here the path of the folder where the output of the Fiji scripts is saved
fijiOutFolder = '/Volumes/Lionnet1/20230812-matlab/20230812-IntS11ivsCtrl_Ser5ph-out';
fijiOutFolder = '/Volumes/lionnt01lab/lionnt01labspace/Feiyue_Tim/20230812-Ser5ph dataset with IntS11i vs control, FP vs control/20230812-IntS11ivsCtrl_Ser5ph-out';

% entere here the path where to re-load the data from a previous Matlab run
% (only used if useMatlabRatherThanFiji = 1;)
matlabInFolder = '/Volumes/lionnt01lab/lionnt01labspace/Feiyue_Tim/20230812-Ser5ph dataset with IntS11i vs control, FP vs control/20230812-IntS11ivsCtrl_Ser5ph-matlabout';
matlabInFolder = '/Users/lionnt01/Dropbox/data/feiyue/Ints11i_Ser5ph_TRI/matlabOut2';

% entere here the path where to save the data 
matlabOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/Ints11i_Ser5ph_TRI/matlabOut2';

% if the data has already been uploaded in Matlab and saved in the
% matlabOutFolder, set this Flag to 1 (faster)
useMatlabRatherThanFiji = 1;

% (Only used if useMatlabRatherThanFiji = 1) 
% If set to 1, the following variable will override the matlabOutFolder folder/subfolders in the
% .mat file, and replace them with folder/subfolders based on the
% matlabOutFolder hardcoded above.
overrideMatFileOutputPaths = 1;

% conditions order
    % conditions are listed in ec.conditionNames. this 
    % array sets the order in which they will be plotted 
conditionsOrder = [1,2,3,4]; 

% set path for subfunctions
addpath('subfunctions/');
addpath('cbrewer/');

% channel names for displays
channelNamesForDisplays = {'DAPI','Pol2','Ser5ph','P-TEFb'};

% eggChamberStages to Include in plots
ecStagesToInclude = 0:10;

% HLB min/max volume in um^3
hlbMinVol = 1;
hlbMaxVol = Inf;

% small cluster min/max volume
smallClusterMinVol = hlbMinVol/50;
smallClusterMaxVol = hlbMinVol/5;

% limitPlots
limitPlots = 1;
%% upload the data

if useMatlabRatherThanFiji
    disp('loading data from matlab saved workspace...');
    if overrideMatFileOutputPaths
        hardCodedMatlabOutFolder = matlabOutFolder;
    end
    load(fullfile(matlabInFolder,'globalAnalysis.mat'));
    if overrideMatFileOutputPaths
        matlabOutFolder = hardCodedMatlabOutFolder;
    end
    disp('done.');
else

    % intitialize eggChamberDataFolder objec - this collecs the locations of
        % the files associated with the experiment and stores them by condition and sample
    ec = eggChamberDataSet(fijiOutFolder);
    
   % load nucleus-wide raw data
    % load raw data from nuclei segmentation - builds a table e
        % c.nucFullT that combines
        % nuclei from all samples (one row per nucleus) and compiles dozens of
        % metrics of size and intensity in every channel. (creates nucT as a backup copy)
    ec.loadAllEggChamberNucleiData;
    
    % load cluster data
        % populates the table ec.clustT
        % also adds nucleoli and nucleoplasm data to ec.nucT and ec.nucFullT tables
        % does NOT add summary cluster metrics to nucT and nucFullT
    ec.loadAllClusterData();
    
    % background correct nuclei intensity
        % subtracts image-wide background value (wholeImgSubtracted) and typical background
        % in an eggchamber away from nuclei (sampleROISubtracted)
    ec.backgroundCorrectNucIntensity;
  
    % add nuc stats to cluster table
        % takes metrics from nucFullT and copies them as extra variables of the
        % cluster table. 
        
    disp('adding Nuclei stats to Cluster table...');  
    ec.clustT = ec.addNucStatsToClustTable();
    
    % background correct clusters intensity
        % subtracts nucleoplasm intensity (plasmSubtracted) and nucleoli intensity
        % (nucleoliSubtracted)
    ec.backgroundCorrectClustIntensity();

    % add cluster data to nuc table
        % takes metrics from clustT, averages them over each nucleus
        % and copies them as extra variables of the nuc table.
    disp('adding cluster stats to nuclei table...');    
    % small clusters
    ec.addAverageClusterStatsToNucTable(smallClusterMinVol,smallClusterMaxVol);

    % large clusters
    ec.addAverageClusterStatsToNucTable(hlbMinVol,Inf);  
    disp('done');
end

%% save workspace in Matlab Output Folder, prepare subfolders
if ~exist(matlabOutFolder,'dir')
    mkdir(matlabOutFolder);
end
save(fullfile(matlabOutFolder,'globalAnalysis.mat')); 

nucFolder = fullfile(matlabOutFolder,'nucleiMetrics');
if ~exist(nucFolder,'dir')
    mkdir(nucFolder);
end

clustFolder = fullfile(matlabOutFolder,'clusterMetrics');
if ~exist(clustFolder,'dir')
    mkdir(clustFolder);
end

qcFolder = fullfile(matlabOutFolder,'qcMetrics');
if ~exist(qcFolder,'dir')
    mkdir(qcFolder);
end

%% set up some plotting items

% colors for the conditions
condColors = cbrewer('qual','Dark2',numel(ec.condIndices));
if numel(ec.condIndices) == 2
    condColors = condColors([1,3],:);
end

% exclude nuclei and clusters outside of the eggChamber of interest
nucleiInGoodEggChambers = ec.nucFullT.eggChamber_Label ~=0 ;
clustersInGoodEggChambers = ec.clustT.eggChamber_Label ~=0 ;

% build name to add to metrics for small and large clusters
smallClusterMetricName = [num2str(smallClusterMinVol),'_',num2str(smallClusterMaxVol)];
largeClusterMetricName = [num2str(hlbMinVol),'_',num2str(hlbMaxVol)];

%% scatter plot by egg chamber: Nucleus Volume (qc)
% plot where all nuclei from an egg chamber are grouped separately.

% all stages
[fullTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
    'nuc',1,'Volume','','nucVolume',nucleiInGoodEggChambers,{},...
    conditionsOrder,ecStagesToInclude,1,'useMean',0.1);

saveNucDataFromPlot(fh,fullTable,avgEcTable,avgCondTable, qcFolder,'nucVolume');

%% scatter plot nuclei intensity by egg chamber (qc): 
% these are good for QC pourposes but not the best metric - use the
% nucleoplasm intensity computed in the next code block instead.

% raw:
for i=1:numel(channelNamesForDisplays)
    [nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
        'nuc',i,'Mean','raw',...
        ['meanNucInt',channelNamesForDisplays{i},'_raw'],nucleiInGoodEggChambers,...
        {},conditionsOrder,ecStagesToInclude,1,'useMean',0.3);

    saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, qcFolder,...
        ['meanNucInt',channelNamesForDisplays{i},'_raw']);
end

% eggChamber-corrected:
for i=1:numel(channelNamesForDisplays)

    [nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
        'nuc',i,'Mean','eggChamberSubtracted',...
        ['meanNucInt',channelNamesForDisplays{i},'_eggChamberSubtracted'],nucleiInGoodEggChambers,...
        {},conditionsOrder,ecStagesToInclude,1,'useMean',0.3);

    saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, qcFolder,...
        ['meanNucInt',channelNamesForDisplays{i},'_eggChamberSubtracted']);
end

%% scatter plot nucleoplasm intensity by egg chamber (nucleoli subtracted): 
% this is the best metric for nuclear expression levels

% nucleoplasm corrected for background estimated across the nucleolus:
for i=1:numel(channelNamesForDisplays)

    [nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
        'plasm',i,'Mean','nucleoliSubtracted',...
        ['meanPlasmInt',channelNamesForDisplays{i},'_nucleoliSub'],nucleiInGoodEggChambers,...
        {},conditionsOrder,ecStagesToInclude,1,'useMean',0.3);

    saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, nucFolder,...
        ['meanPlasmInt',channelNamesForDisplays{i},'_nucleoliSub']);
end

%% cluster intensity as a function of cluster volume, raw (qc)
scatterPlotClusterIntVsCluterVolume(ec,ecStagesToInclude,clustersInGoodEggChambers,...
    conditionsOrder,condColors,channelNamesForDisplays,'raw',qcFolder);

%% cluster intensity as a function of cluster volume, normalized by nucleoplasm intensity (qc)
scatterPlotClusterIntVsCluterVolume(ec,ecStagesToInclude,clustersInGoodEggChambers,...
    conditionsOrder,condColors,channelNamesForDisplays,'plasmNorm',qcFolder);

%% scatter plot by egg chamber: Number of large clusters per nucleus, Cluster volume

nLargeClustersMetricName = strrep(['NumClusters',largeClusterMetricName],'.','pt');

% plot number of large clusters per nucleus
[nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
    'nuc',1,nLargeClustersMetricName,'','nLargeClusters',nucleiInGoodEggChambers,...
    {},conditionsOrder,ecStagesToInclude,1,'useMean',0.2);
%
saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, nucFolder,'nLargeClusters');

% plot cluster volume for large clusters
[clustTable,avgNucClustTable,avgEcTable,avgCondTable,fhFull] = ec.scatterPlotClustTableMetricByEggChamber( ...
    'clust',1,'Volume','','clusterVolume',clustersInGoodEggChambers,...
    {},conditionsOrder,ecStagesToInclude,hlbMinVol,hlbMaxVol,1,'useMean',0.1);

% same, but nucleus average
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'mean_clusterVolume',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'mean large cluster Volume by Nuclei',0.1);

% save tables
saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcTable,avgCondTable, clustFolder,'clustVolume');


%% plot number of cluster per nucleus vs cluster intensity, per channel
doPlot = 1; % set to 0 to skip
if doPlot
    idx0 = ismember(ec.nucFullT.eggChamber_Stage , ecStagesToInclude) ...
            & nucleiInGoodEggChambers;
        
    for i=1:numel(channelNamesForDisplays)
        figure('Name',['nClusters/nucleus vs cluster intensity ',channelNamesForDisplays{i}]);
        hold on;
        for j=1:numel(ec.condIndices)
                curColor = condColors(j,:);
                curIdx = idx0 & ec.nucFullT.cond_Idx == conditionsOrder(j);
                curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
                xData = ec.nucFullT.(['nuc_NumClusters',largeClusterMetricName]); % number of cluster in the nucleus
                yData = ec.nucFullT.(['nucAvgClust',largeClusterMetricName,'_clust_C',num2str(i),'Median_nucleoliSubtracted'])...
                    ./ec.nucFullT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);
                yData(xData == 0) = 2.^(log(max(yData(~isnan(yData)))-min(yData(~isnan(yData))))/log(2)*(rand(sum(xData ==0),1)));
                yData = log(yData)/log(10);
                s = scatter(xData(curIdx),yData(curIdx),100,...
                            curColor,'o','filled','DisplayName',ec.conditionNames{conditionsOrder(j)});
                alpha(s,0.1); 
        end
        xlabel('Number of clusters per nucleus');
        ylabel('average cluster intensity, plasmNorm');
    end
end

%% plot number of cluster per um^3 vs cluster intensity, per channel
doPlot = 1; % set to 0 to skip
if doPlot
    idx0 = ismember(ec.nucFullT.eggChamber_Stage , ecStagesToInclude) ...
        & nucleiInGoodEggChambers;
    
    for i=1:numel(channelNamesForDisplays)
        figure('Name',['nClusters/nucleus vs cluster intensity ',channelNamesForDisplays{i}]);
        hold on;
        for j=1:numel(ec.condIndices)
                curColor = condColors(j,:);
                curIdx = idx0 & ec.nucFullT.cond_Idx == conditionsOrder(j);
                curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
                xData = ec.nucFullT.(['nuc_NumClusters',largeClusterMetricName])...
                    ./ec.nucFullT.nuc_Volume; % number of cluster in the nucleus divided by nuc volume
                yData = ec.nucFullT.(['nucAvgClust',largeClusterMetricName,'_clust_C',num2str(i),'Median_nucleoliSubtracted'])...
                    ./ec.nucFullT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);
                yData(xData == 0) = 2.^(log(max(yData(~isnan(yData)))-min(yData(~isnan(yData))))/log(2)*(rand(sum(xData ==0),1)));
                yData = log(yData)/log(2);
                s = scatter(xData(curIdx),yData(curIdx),100,...
                            curColor,'o','filled','DisplayName',ec.conditionNames{conditionsOrder(j)});
                alpha(s,0.1);
        end
        xlabel('Number of clusters per unit volume');
        ylabel('average cluster intensity, plasmNorm');
    end
end
%% cdf of cluster intensity by nucleus
cdfBinSize = 0.02;
idx0 = ismember(ec.nucFullT.eggChamber_Stage , ecStagesToInclude) ...
        & nucleiInGoodEggChambers;
    
for i=1:numel(channelNamesForDisplays)
    fh = figure('Name',['cdf: cluster intensity ',channelNamesForDisplays{i},' plasmNorm by nucleus']);
    hold on;
    yData = ec.nucFullT.(...
        ['nucAvgClust',largeClusterMetricName,'_clust_C',num2str(i),'Median_nucleoliSubtracted'])...
        ./ec.nucFullT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);
    nInf = 0;
    nNeg = 0;
    nTot = 0;
    % find min and max values to plot all cdfs along the same x values 
    % (makes it easier to save as single table)
    jIdx = idx0 & ~isnan(yData) & ~isinf(yData) & yData>0;
    condMin = min(log(yData(jIdx))/log(10));
    condMax = max(log(yData(jIdx))/log(10));
    
    for j=1:numel(ec.condIndices)
        curColor = condColors(j,:);
        curIdx = idx0 & ec.nucFullT.cond_Idx == conditionsOrder(j);
        curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
        
        curData = yData(curIdx);
        nTot = nTot+numel(curData);

        % starting number of datapoints
        n0 = numel(curData);

        % get NaN statistics (these are the nuclei without clusters)
        percentEmptyNuclei = sum(isnan(curData))/numel(curData)*100;
        nEmptyNuclei = sum(isnan(curData));
        
        curData = curData(~isnan(curData));
        n1 = numel(curData);
        % remove negative values before taking the log
        curData = curData(curData>0);
        n2 = numel(curData);

        % removing any inf values 
        curData = curData(~isinf(curData));
        n3 = numel(curData);

        % generating stats for warning msg in case some negative/inf values were
        % removed
        nInf = nInf+ n2-n3;
        nNeg = nNeg + n1-n2;

        % taking the log of the cleaned up dataset.
        curData = log(curData)/log(10);
        
        [n,x] = hist(curData,condMin:cdfBinSize:condMax);
        n(1) = n(1) + nEmptyNuclei; % adding empty nuclei to first bin
        plot(x,cumsum(n)/sum(n),...
            'Color',curColor,'LineWidth',2,'DisplayName',...
            [ec.conditionNames{conditionsOrder(j)},', ',num2str(percentEmptyNuclei,2),'% nuclei without detected clusters']);

        % add data to table to be saved
        if j==1
            curT = table(x',cumsum(n')/sum(n),'VariableNames',...
                {['log10nucAvgClust_C',num2str(i),'Median_nucleoliSub_plasmNorm'],...
                [ec.conditionNames{conditionsOrder(j)},'_Count']});
        else
            curT =addvars(curT,cumsum(n')/sum(n),'NewVariableNames',{[ec.conditionNames{conditionsOrder(j)},'_Count']});
        end
    end
    warningMsg = {};
    if nNeg > 0
        warningMsg{end+1} = ['Warning: discarded ',num2str(nNeg),'/',num2str(nTot),' datapoints with negative values.'];
    end
    if nInf > 0
        warningMsg{end+1} = ['Warning: discarded ',num2str(nInf),'/',num2str(nTot),' datapoints with infinite values.'];
    end

    if ~isempty(warningMsg)
        dim = [.5 .1 .3 .3];
        annotation('textbox',dim,'String',warningMsg,'FitBoxToText','on',...
            'Units','normalized','BackgroundColor','white');
    end
    xlabel(['C',num2str(i),' avg cluster Int, plasm norm, log10']);
    ylabel('Cumulated Distribution Function');
    grid on;
    legend show;

    % save figure and data
    figName = ['cdf_nucAvgClust_C',num2str(i),'Median_nucleoliSub_plasmNorm'];
    saveas(fh,fullfile(clustFolder,[figName,'.fig']));
    try
        saveas(fh,fullfile(clustFolder,[figName,'.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',figName,' - likely permission issue.']);
    end
    writetable(curT,fullfile(clustFolder,[figName,'.txt']),'Delimiter','\t');
end

%% plot distribution of the number of HLBs per nucleus in G vs S
GvsSChannel = 4;
GvsSThresh = [0.45, 0.45, 0.66, 0.66]; % per condition, threshold based on the previous section curves (log10 values)
histBinSize = 1;
varName = ['nuc_NumClusters',largeClusterMetricName];
chIdx = 0;
rawOrPlasmNorm = 'raw';
logTransform = 0;
histXLabel = 'Number of HLBs detected per nucleus';
filePrefix = 'HistGvsS';

plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
    nucleiInGoodEggChambers,conditionsOrder,condColors,histBinSize,...
    channelNamesForDisplays,chIdx,rawOrPlasmNorm,nucFolder,...
    logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName);

%% plot histogram of the number of small clusters per nucleus in G vs S
GvsSChannel = 4;
GvsSThresh = [0.45, 0.45, 0.66, 0.66]; % per condition, threshold based on the previous section curves (log10 values)
histBinSize = 10;
varName = ['nuc_NumClusters',strrep(smallClusterMetricName,'.','pt')];
chIdx = 0;
rawOrPlasmNorm = 'raw';
logTransform = 0;
histXLabel = 'Number of small Clusters detected per nucleus';
filePrefix = 'HistGvsS';

plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
    nucleiInGoodEggChambers,conditionsOrder,condColors,histBinSize,...
    channelNamesForDisplays,chIdx,rawOrPlasmNorm,nucFolder,...
    logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName);

%% scatter plot by egg chamber: Number of small clusters per nucleus, Cluster volume

% plot cluster volume for small clusters
[clustTable,avgNucClustTable,avgEcTable,avgCondTable,fhFull] = ec.scatterPlotClustTableMetricByEggChamber( ...
    'clust',1,'Volume','','clusterVolume',clustersInGoodEggChambers,...
    {},conditionsOrder,ecStagesToInclude,smallClusterMinVol,smallClusterMaxVol,1,'useMean',0.01);

% same, but nucleus average
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'mean_clusterVolume',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'mean small cluster Volume by Nuclei',0.05);

% plot the number of small clusters per nucleus
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'nClusters',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'number of Small Clusters per Nucleus',0.05);

% save tables
saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcTable,avgCondTable, clustFolder,'smallClustVolume');

%% cluster intensity, HLBs (raw)
if ~limitPlots
    for i=1:numel(channelNamesForDisplays)
        curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
        [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
            ec.scatterPlotClustTableMetricByEggChamber(...
            'clust',i,'Median','nucleoliSubtracted',curVarName,...
            clustersInGoodEggChambers,{},conditionsOrder,ecStagesToInclude,...
            hlbMinVol,hlbMaxVol,1,'useMean',0.05);
    
        fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
            ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
            conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);
    
        saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
            ['clustInt',channelNamesForDisplays{i},'_nucleoliSub']);
    end
end
%% histogram of HLB cluster volumes by condition
histBinSize = 0.03;
logTransform = 1;
rawOrPlasmNorm = 'raw';
histXLabel = 'Cluster Volume in um^3';
filePrefix = 'HLB';
minVol = hlbMinVol; maxVol = hlbMaxVol;

plotClusterVarHistogram(ec,'clust_Volume',ecStagesToInclude,...
    clustersInGoodEggChambers,conditionsOrder,condColors,...
    minVol,maxVol,histBinSize,...
    channelNamesForDisplays,0,rawOrPlasmNorm,clustFolder,logTransform,histXLabel,filePrefix);

%% plot histogram of the HLB volume in G vs S
GvsSChannel = 4;
GvsSThresh = [0.45, 0.45, 0.66, 0.66]; % per condition, threshold based on the previous section curves (log10 values)
histBinSize = 1;
varName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_Volume'];
chIdx = 0;
rawOrPlasmNorm = 'raw';
logTransform = 0;
histXLabel = 'Average Cluster Volume per nucleus (um^3)';
filePrefix = 'HistGvsS';

plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
    nucleiInGoodEggChambers,conditionsOrder,condColors,histBinSize,...
    channelNamesForDisplays,chIdx,rawOrPlasmNorm,nucFolder,...
    logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName);

%% plot histograms of the HLB intensity in G vs S

GvsSChannel = 4;
GvsSThresh = [0.45, 0.45, 0.66, 0.66]; % per condition, threshold based on the previous section curves (log10 values)
histBinSize = 0.1;
rawOrPlasmNorm = 'plasmNorm';
logTransform = 1;
if logTransform
    logLabel = ' log10';
else
    logLabel = '';
end
filePrefix = 'HistGvsS';

for chIdx=1:numel(channelNamesForDisplays)
    varName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_C',num2str(chIdx),'Median_nucleoliSubtracted'];
    histXLabel = ['Average HLB ',channelNamesForDisplays{chIdx},' Intensity nucleoliSub ',rawOrPlasmNorm,logLabel];

    plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
        nucleiInGoodEggChambers,conditionsOrder,condColors,histBinSize,...
        channelNamesForDisplays,chIdx,rawOrPlasmNorm,nucFolder,...
        logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName);
end

%% load images of nuclei
imageJViz = 1;

metricChIdx = 4;
metricRange = 10.^[0.2,0.3];
metricRawOrPlasmNorm = 'plasmNorm';
removeBackgroundEggChambers = 1;

% list nuclei rows that fit in a range of a nucleus metric
%idxData = listNucleiWithinRange(metricChIdx,metricRange,metricRawOrPlasmNorm);
switch metricRawOrPlasmNorm
    case 'raw'
        metricVarName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_C',num2str(metricChIdx),'Median_nucleoliSubtracted'];
        yData = ec.nucFullT.(metricVarName);
    case 'plasmNorm'
        metricVarName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_C',num2str(metricChIdx),'Median_nucleoliSubtracted'];
        denomVarName = ['plasm_C',num2str(metricChIdx),'Median_nucleoliSubtracted'];
        yData = ec.nucFullT.(metricVarName)./ec.nucFullT.(denomVarName);
end
idxData = yData >= metricRange(1) & yData <= metricRange(2);
if removeBackgroundEggChambers
    idxData = idxData & ec.nucFullT.eggChamber_Idx ~=0;
end

% generate basic variables
varList = {'sample_InputFileName','cond_Idx','sample_Idx','eggChamber_Idx',...
    'eggChamber_Stage','nuc_Label','nuc_Volume',...
    ['nuc_NumClusters',strrep(largeClusterMetricName,'.','pt')],...
    ['nuc_NumClusters',strrep(smallClusterMetricName,'.','pt')]};

varIdx = ismember(ec.nucFullT.Properties.VariableNames,varList);

nucIdxList = ec.nucFullT(idxData,varIdx);

% add channel main metrics
for i=1:numel(channelNamesForDisplays)
    curVarName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_C',num2str(i),'Median_raw'];
    nucIdxList = addvars(nucIdxList,ec.nucFullT.(curVarName)(idxData),'NewVariableNames',{curVarName});
    
    curVarName = ['plasm_C',num2str(i),'Median_raw'];
    nucIdxList = addvars(nucIdxList,ec.nucFullT.(curVarName)(idxData),'NewVariableNames',{curVarName});
    
    curVarName = ['nucleoli_C',num2str(i),'Median_raw'];
    nucIdxList = addvars(nucIdxList,ec.nucFullT.(curVarName)(idxData),'NewVariableNames',{curVarName});

    metricVarName = ['nucAvgClust',strrep(largeClusterMetricName,'.','pt'),'_clust_C',num2str(i),'Median_nucleoliSubtracted'];
    denomVarName = ['plasm_C',num2str(i),'Median_nucleoliSubtracted'];
    yData = ec.nucFullT.(metricVarName)./ec.nucFullT.(denomVarName);
    nucIdxList = addvars(nucIdxList,yData(idxData),'NewVariableNames',{[metricVarName,'_plasmCorr']});

end

% load n example nuclei
if imageJviz
    nNucleiToLoad = 5;
    idxToLoad = ceil(size(nucIdxList,1)*rand(nNucleiToLoad,1));
    % load IJ interfacte
    MIJ.start;
    
    for i=1:numel(idxToLoad)
        '182-128.tif'
        % reconstruct file name
        nucFileName = fullfile(fijiOutFolder,);

        % loadImage
        image = mijread(nucFileName);
    end

end


%%%%%%%%%%%%%%%%%%%%%
%% plot histograms of small clusters intensity in G vs S

GvsSChannel = 4;
GvsSThresh = [0.45, 0.45, 0.66, 0.66]; % per condition, threshold based on the previous section curves (log10 values)
histBinSize = 0.1;
rawOrPlasmNorm = 'plasmNorm';
logTransform = 1;
if logTransform
    logLabel = ' log10';
else
    logLabel = '';
end
filePrefix = 'HistGvsS';

for chIdx=1:numel(channelNamesForDisplays)
    varName = ['nucAvgClust',strrep(smallClusterMetricName,'.','pt'),'_clust_C',num2str(chIdx),'Median_nucleoliSubtracted'];
    histXLabel = ['Average Small Cluster ',channelNamesForDisplays{chIdx},' Intensity nucleoliSub ',rawOrPlasmNorm,logLabel];

    plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
        nucleiInGoodEggChambers,conditionsOrder,condColors,histBinSize,...
        channelNamesForDisplays,chIdx,rawOrPlasmNorm,nucFolder,...
        logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName);
end

%% scatter plot of cluster intensity for SMALL clusters (nucleoli subtracted)
if ~limitPlots
    for i=1:numel(channelNamesForDisplays)
       curVarName = ['SmallClustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
        [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
            ec.scatterPlotClustTableMetricByEggChamber(...
            'clust',i,'Median','nucleoliSubtracted',curVarName,...
            clustersInGoodEggChambers,{},conditionsOrder,ecStagesToInclude,...
            smallClusterMinVol,smallClusterMaxVol,1,'useMean',0.01);
    
        fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
            ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
            conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.05);
    
        saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
            ['smallClustInt',channelNamesForDisplays{i},'_nucleoliSub']);
    end
end
%% scatter plot cluster intensity, HLBs (nucleoli subtracted, relative to plasm levels)

for i=1:numel(channelNamesForDisplays)
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);
    
    curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtr_plasmNorm'];

    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustArbitraryMetricByEggChamber(yData,curVarName,...
        clustersInGoodEggChambers,{},conditionsOrder,ecStagesToInclude,...
        hlbMinVol,hlbMaxVol,1,'useMean',0.2);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['clustInt_',channelNamesForDisplays{i},'_plasmNorm']);
end

%% scatter plot of cluster intensity for SMALL clusters (nucleoli subtracted, relative to plasm levels)
for i=1:numel(channelNamesForDisplays)
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);

    curVarName = ['smallClustInt_',channelNamesForDisplays{i},'_nucleoliSubtr_plasmNorm'];

    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustArbitraryMetricByEggChamber(...
        yData,curVarName,...
        true(size(ec.clustT,1),1),{},conditionsOrder,ecStagesToInclude,...
        smallClusterMinVol,smallClusterMaxVol,1,'useMean',0.01);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['smallClustInt_',channelNamesForDisplays{i},'_plasmNorm']);
end

%% histogram of HLB cluster intensities (relative to plasm levels) by condition 

    histBinSize = 0.05;
    logTransform = 1;
    rawOrPlasmNorm = 'plasmNorm';
    filePrefix = 'HLB';
    minVol = hlbMinVol; maxVol = hlbMaxVol;
    
    for i=1:numel(channelNamesForDisplays)
        histXLabel = ['HLB C',num2str(i),' Median nucleoliSub, ',rawOrPlasmNorm];
        
        plotClusterVarHistogram(ec,...
            ['clust_C',num2str(i),'Median_nucleoliSubtracted'],ecStagesToInclude,...
            clustersInGoodEggChambers,conditionsOrder,condColors,...
            minVol,maxVol,histBinSize,channelNamesForDisplays,i,rawOrPlasmNorm,...
            clustFolder,logTransform,histXLabel,filePrefix);
    end

%% histogram of SMALL cluster intensities (relative to plasm levels) by condition
if ~limitPlots
    histBinSize = 0.05;
    logTransform = 0;
    rawOrPlasmNorm = 'plasmNorm';
    filePrefix = 'smallClust';
    minVol = smallClusterMinVol; maxVol = smallClusterMaxVol;
    
    for i=1:numel(channelNamesForDisplays)
        histXLabel = ['small clust C',num2str(i),' Median nucleoliSub, ',rawOrPlasmNorm];
        
        plotClusterVarHistogram(ec,...
            ['clust_C',num2str(i),'Median_nucleoliSubtracted'],ecStagesToInclude,...
            clustersInGoodEggChambers,conditionsOrder,condColors,...
            minVol,maxVol,histBinSize,channelNamesForDisplays,i,rawOrPlasmNorm,...
            clustFolder,logTransform,histXLabel,filePrefix);
    end
end
%% histogram of HLB cluster intensities (raw) by condition
if ~limitPlots
    histBinSize = 100;
    logTransform = 0;
    rawOrPlasmNorm = 'raw';
    filePrefix = 'HLBClust';
    minVol = hlbMinVol; maxVol = hlbMaxVol;
    
    for i=1:numel(channelNamesForDisplays)
        
        histXLabel = ['HLB clust C',num2str(i),' Median nucleoliSub, ',rawOrPlasmNorm];
        
        plotClusterVarHistogram(ec,...
            ['clust_C',num2str(i),'Median_nucleoliSubtracted'],ecStagesToInclude,...
            clustersInGoodEggChambers,conditionsOrder,condColors,...
            minVol,maxVol,histBinSize,channelNamesForDisplays,i,rawOrPlasmNorm,...
            clustFolder,logTransform,histXLabel,filePrefix);
    end
end
%% histogram of SMALL cluster intensities (raw) by condition
if ~limitPlots
    histBinSize = 50;
    logTransform = 0;
    rawOrPlasmNorm = 'raw';
    filePrefix = 'HLBClust';
    minVol = smallClusterMinVol; maxVol = smallClusterMaxVol;
    
    for i=1:numel(channelNamesForDisplays)
        
        histXLabel = ['HLB clust C',num2str(i),' Median nucleoliSub, ',rawOrPlasmNorm];
        
        plotClusterVarHistogram(ec,...
            ['clust_C',num2str(i),'Median_nucleoliSubtracted'],ecStagesToInclude,...
            clustersInGoodEggChambers,conditionsOrder,condColors,...
            minVol,maxVol,histBinSize,channelNamesForDisplays,i,rawOrPlasmNorm,...
            clustFolder,logTransform,histXLabel,filePrefix);  
    end
end
%% HLB intensity, channel vs channel (nucleoli subtracted, plasm Norm)
chX = 4;  % channel for X axis
chY = 2;  % channel for Y axis
alphaVal = 0.3; % transparency
rawOrPlasmNorm = 'plasmNorm';
varPrefix = 'clust';
varBaseName = 'Median';
varSuffix = 'nucleoliSubtracted';
filePrefix = 'hlb';
scatterClusterVarChannelVsChannel(ec,varPrefix,varBaseName,varSuffix,...
    rawOrPlasmNorm,clustersInGoodEggChambers,chX,chY,hlbMinVol,hlbMaxVol,...
    ecStagesToInclude,conditionsOrder,condColors,channelNamesForDisplays,alphaVal,qcFolder,filePrefix);

%% Small Cluster intensity, channel vs channel (nucleoli subtracted, plasm Norm)
chX = 4; % channel for X axis
chY = 2; % channel for Y axis
alphaVal = 0.01; % transparency
rawOrPlasmNorm = 'plasmNorm';
varPrefix = 'clust';
varBaseName = 'Median';
varSuffix = 'nucleoliSubtracted';
filePrefix = 'smallClust';
scatterClusterVarChannelVsChannel(ec,varPrefix,varBaseName,varSuffix,...
    rawOrPlasmNorm,clustersInGoodEggChambers,chX,chY,smallClusterMinVol,smallClusterMaxVol,...
    ecStagesToInclude,conditionsOrder,condColors,channelNamesForDisplays,alphaVal,qcFolder,filePrefix);

%% HLB chY/chX ratio, by sample - Normalized to nuclear levels (nucleoli Subtracted, normalized to nucleoplasm levels)
chY = 3; % numerator channel
chX = 2; % denominator channel

minVolume = hlbMinVol; % minimum cluster Volume
maxVolume = hlbMaxVol; % max cluster Volume

% egg Chambers/clusters to include based on filters
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0 ...
    & ec.clustT.clust_Volume >= minVolume ...
    & ec.clustT.clust_Volume <= maxVolume;

numData = ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted']);

denomData = ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted']); 

yData = numData./denomData;

[clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fh] = ...
    ec.scatterPlotClustArbitraryMetricByEggChamber(...
    yData,[channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX},', normalized by nucleoplasm'],idx0,...
    {},conditionsOrder,ecStagesToInclude,minVolume,maxVolume,1,'useMedian',0.1);

saveNucDataFromPlot(fh,clustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm']);

%% histogram of HLB cluster chY/chX ratio (relative to plasm levels) by condition
chY = 3; % numerator channel
chX = 2; % denominator channel
histBinSize = 0.2;
logTransform = 0;
figTitle = ['Hist HLB ',channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX},'normalized by nucleoplasm'];
histXLabel = ['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'];
filePrefix = '';
varName = ['hlbInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'];

numData = ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted']);
denomData = ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted']); 
yData = numData./denomData;

plotArbitraryClusterVarHistogram(ec,yData,varName,figTitle,ecStagesToInclude,...
    clustersInGoodEggChambers,conditionsOrder,condColors, hlbMinVol,hlbMaxVol,histBinSize,...
    clustFolder,logTransform,histXLabel,filePrefix);

%% histogram of SMALL cluster chY/chX ratio (relative to plasm levels) by condition
chY = 3; % numerator channel
chX = 2; % denominator channel
histBinSize = 0.2;
logTransform = 0;
figTitle = ['Hist small Clusters ',channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX},'normalized by nucleoplasm'];
histXLabel = ['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'];
filePrefix = '';
varName = ['smallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'];

numData = ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted']);
denomData = ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted']); 
yData = numData./denomData;

plotArbitraryClusterVarHistogram(ec,yData,varName,figTitle,ecStagesToInclude,...
    clustersInGoodEggChambers,conditionsOrder,condColors,smallClusterMinVol,smallClusterMaxVol,histBinSize,...
    clustFolder,logTransform,histXLabel,filePrefix);

%% plot cluster intensity, stratified by ptefb

refPrefix = 'clust';
refChannel = 4;
refBaseName = 'Median';
refSuffix = 'nucleoliSubtracted';
refVarType = 'channel';
refRawOrPlasmNorm = 'plasmNorm';
dataPrefix = 'clust';
dataBaseName = 'Median';
dataSuffix = 'nucleoliSubtracted';
dataVarType = 'channel';
dataRawOrPlasmNorm = 'plasmNorm';
yMin0 = [-Inf,-Inf,-Inf,-Inf];
yMax0 = [3,3,3,4.5];
yMin1 = yMax0;
yMax1= [Inf,Inf,Inf,Inf];
stratificationName = 'PTEFBpositive';
fileSuffix = 'ptefbStratified';
alphaVal = 0.3;
for i=1:numel(channelNamesForDisplays)
    dataChannel = i;
    plotClusterVarStratified(ec,refPrefix,refChannel,refBaseName,refSuffix,refVarType,refRawOrPlasmNorm,...
        dataPrefix,dataChannel,dataBaseName,dataSuffix,dataVarType,dataRawOrPlasmNorm,...
        clustersInGoodEggChambers,yMin0,yMax0,yMin1,yMax1,stratificationName,hlbMinVol,hlbMaxVol,fileSuffix,...
        alphaVal,condColors,conditionsOrder,ecStagesToInclude,clustFolder);
end

%% plotting ser5/pol2 ratio, stratified by ptefb
refChannel = 4; % channel used to stratify (ptefb)
chY = 3; % numerator
chX = 2; % denominator

yMin0 = [-Inf,-Inf,-Inf,-Inf];
yMax0 = [3,3,3,4.5];
yMin1 = yMax0;
yMax1= [Inf,Inf,Inf,Inf];
stratificationName = 'PTEFBpositive';
fileSuffix = 'ptefbStratified';
alphaVal = 0.3;

yRef = ec.clustT.(['clust_C',num2str(refChannel),'Median_nucleoliSubtracted']) ...
    ./ec.clustT.(['plasm_C',num2str(refChannel),'Median_nucleoliSubtracted']);

numData = ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted']);

denomData = ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted']); 

yData = numData./denomData;

plotClusterArbitraryVarStratified(ec,yRef,yData,clustersInGoodEggChambers,yMin0,yMax0,yMin1,yMax1,...
    stratificationName,'ser5ph_pol2_ratio',hlbMinVol,hlbMaxVol,fileSuffix,alphaVal,condColors,...
    conditionsOrder,ecStagesToInclude,clustFolder);

%% plotting MPM2+ and MPM2- HLB Volume side by side

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting functions

%%
function plotClusterVarStratified(ec,refPrefix,refChannel,refBaseName,refSuffix,refVarType,refRawOrPlasmNorm,...
    dataPrefix,dataChannel,dataBaseName,dataSuffix,dataVarType,dataRawOrPlasmNorm,...
    idxData,yMin0,yMax0,yMin1,yMax1,stratificationName,minVolume,maxVolume,fileSuffix,...
    alphaVal,condColors,conditionsOrder,ecStagesToInclude,clustFolder)
    
    refVar = ec.buildVarName(refPrefix,refChannel,refBaseName,refSuffix,refVarType);
    switch refRawOrPlasmNorm
        case 'raw'
            yRef = ec.clustT.(refVar);
        case 'plasmNorm'
            refVarDenom = ec.buildVarName('plasm',refChannel,refBaseName,refSuffix,refVarType);
            yRef = ec.clustT.(refVar)./ec.clustT.(refVarDenom);
    end
    
    dataVar = ec.buildVarName(dataPrefix,dataChannel,dataBaseName,dataSuffix,dataVarType);
    switch dataRawOrPlasmNorm
        case 'raw'
            yData = ec.clustT.(dataVar);
        case 'plasmNorm'
            dataVarDenom = ec.buildVarName('plasm',dataChannel,dataBaseName,dataSuffix,dataVarType);
            yData = ec.clustT.(dataVar)./ec.clustT.(dataVarDenom);
    end

    plotClusterArbitraryVarStratified(ec,yRef,yData,idxData,yMin0,yMax0,yMin1,yMax1,...
        stratificationName,dataVar,minVolume,maxVolume,fileSuffix,alphaVal,condColors,...
        conditionsOrder,ecStagesToInclude,clustFolder);
end
%% plot cluster arbitrary variable, stratified
function plotClusterArbitraryVarStratified(ec,yRef,yData,idxData,yMin0,yMax0,yMin1,yMax1,...
    stratificationName,dataVarName,minVolume,maxVolume,fileSuffix,alphaVal,condColors,...
    conditionsOrder,ecStagesToInclude,clustFolder)

    % chRef: channel used to startify the data in positive vs negative
    % yMin0  min value to be called a negative per condition
    % yMax0  max value to be called a negative per condition
    % yMin1  min value to be called a negative per condition
    % yMax1  max value to be called a negative per condition

    idx0StratificationVal = 0; % negative
    idx1StratificationVal = 1; % positive
    
    if numel(yMin0) ~= ec.nConditions || numel(yMax0) ~= ec.nConditions || numel(yMin1) ~= ec.nConditions || numel(yMax0) ~= ec.nConditions
        disp('yMin0, yMax0, yMin1 and yMax1 need to have the same number of elements as there are conditions');
        return
    end
    idxData0 = false(size(ec.clustT,1),1);
    idxData1 = idxData0;
    for i=1:ec.nConditions
        idxData0 = idxData0 | ...
                idxData & ec.clustT.cond_Idx == i & yRef > yMin0(i) & yRef < yMax0(i);
        idxData1 = idxData1 | ...
                idxData & ec.clustT.cond_Idx == i & yRef > yMin1(i) & yRef < yMax1(i);
    end
    
    % first category (idxData1) is gray, second (idxData2) is in the default colors
    dualColorMap = {0.7*ones(size(condColors));condColors};
    
    [clustTable0,avgNucClustTable0,avgEcClustTable0,avgCondClustTable0,...
        clustTable1,avgNucClustTable1,avgEcClustTable1,avgCondClustTable1,fh] = ...
                ec.scatterPlotClustDualArbitraryMetricByEggChamber(...
                yData,dataVarName,idxData0,idxData1,...
                {},conditionsOrder,ecStagesToInclude,...
                minVolume,maxVolume,1,'useMedian',alphaVal,dualColorMap);
    ylim('auto');

    clustTable0 = addvars(clustTable0,idx0StratificationVal*ones(size(clustTable0,1),1),'newVariableNames',{stratificationName});
    avgNucClustTable0 = addvars(avgNucClustTable0,idx0StratificationVal*ones(size(avgNucClustTable0,1),1),'newVariableNames',{stratificationName});
    avgEcClustTable0 = addvars(avgEcClustTable0,idx0StratificationVal*ones(size(avgEcClustTable0,1),1),'newVariableNames',{stratificationName});
    avgCondClustTable0 = addvars(avgCondClustTable0,idx0StratificationVal*ones(size(avgCondClustTable0,1),1),'newVariableNames',{stratificationName});
    
    clustTable1 = addvars(clustTable1,idx1StratificationVal*ones(size(clustTable1,1),1),'newVariableNames',{stratificationName});
    avgNucClustTable1 = addvars(avgNucClustTable1,idx1StratificationVal*ones(size(avgNucClustTable1,1),1),'newVariableNames',{stratificationName});
    avgEcClustTable1 = addvars(avgEcClustTable1,idx1StratificationVal*ones(size(avgEcClustTable1,1),1),'newVariableNames',{stratificationName});
    avgCondClustTable1 = addvars(avgCondClustTable1,idx1StratificationVal*ones(size(avgCondClustTable1,1),1),'newVariableNames',{stratificationName});

    clustTable = [clustTable0;clustTable1];
    avgNucClustTable = [avgNucClustTable0;avgNucClustTable1];
    avgEcClustTable = [avgEcClustTable0;avgEcClustTable1];
    avgCondClustTable = [avgCondClustTable0;avgCondClustTable1];

    saveClustDataFromPlot(fh,[],clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,[dataVarName,'_',fileSuffix]);
end

%%
function scatterPlotClusterIntVsCluterVolume(ec,ecStagesToInclude,...
    clustersInGoodEggChambers,conditionsOrder,condColors,...
    channelNamesForDisplays,rawOrPlasmNorm,qcFolder)

    
    idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
        & clustersInGoodEggChambers;
    
    for i=1:numel(channelNamesForDisplays)
        switch rawOrPlasmNorm
            case 'raw'
                fh = figure('Name',...
                        ['Cluster Volume vs ',channelNamesForDisplays{i},' intensity (raw)']);
            case 'plasmNorm'
                fh = figure('Name',...
                        ['Cluster Volume vs ',channelNamesForDisplays{i},' intensity (plasmNorm)']);
        end
        hold on;
        for j=1:numel(ec.condIndices)
            curColor = condColors(j,:);
            curIdx = idx0 & ec.clustT.cond_Idx == conditionsOrder(j);
            switch rawOrPlasmNorm
                case 'raw'
                    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']);
                case 'plasmNorm'
                    yData = log(ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted'])...
                         ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']))/log(2);
            end

            s = scatter(ec.clustT.clust_Volume(curIdx),yData(curIdx),16,...
                curColor,'o','filled','DisplayName',ec.conditionNames{conditionsOrder(j)});
            alpha(s,0.2);
        end
        grid on;
        legend show;
        xlabel('Volume (um^3)');
        set(gca,'xscale','log');

        switch rawOrPlasmNorm
                case 'raw'
                    ylabel(['clust C',num2str(i),'Median nucleoliSub']);
                    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'.fig']));
                case 'plasmNorm'
                    ylabel(['log2 clust C',num2str(i),'Median nucleoliSub plasmNorm']);
                    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'.eps']),'epsc');
        end     
    end
end

%% histogram of cluster variable
% plots a histogram of a variable from the clustT table for each condition, all in the same plot.
% also saves the figures and tables in the clustFolder

% INPUT
% ec: eggChamber object
% yData: variable - should be nCluster x 1 vector.
% varName: name of the variable, needs to match a column name of clustT
% figTitle: title of the figure
% ecStagesToInclude: array of stages, e.g. [7,8] or 0:10 to include all
    % eggchambers if the staging wasnt performed.
% idxData: nCluster x 1  logical of which rows to include in the
    % computation
% conditionsOrder: order of plotting the conditions - e.g. [3,4,2,1] will
    % plot condition 3, then 4, then 2 then 1.
% condColors: colors assigned to each condition, Conditions x 3 array
% minVolume/maxVolume: min and maxVolume allowed for the cluster to be included
% hitBinSize: bin size of the histogram
% channelNamesForDisplays: list of channel names in order to use in figure title
% channelIdx: index of the color channel to plot if variable is a channel
    % (will be used in the title of the figure); value will be ignored if 
    % the variable is not a channel variable
% rawOrPlasmNorm: 
    % choose 'raw' to plot the raw cluster variable; 
    % choose 'plasmNorm' to divide the cluster variable by the corresponding 
        % variable computed in the nucleoplasm (by changing the 'clust_' prefix
        % into 'plasm_').
% clustFolder: subfolder where tables/figs will be saved
% logTransform: whether to perform the histogram on a log10 of the data
% histXLabel: label for the x axis
% filePrefix: prefix to insert in the file name, e.g. 'HLB' or
% 'smallClust'
function plotClusterVarHistogram(ec,varName,ecStagesToInclude,...
    idxData,conditionsOrder,condColors,minVolume,maxVolume,histBinSize,...
    channelNamesForDisplays,channelIdx,rawOrPlasmNorm,clustFolder,...
    logTransform,histXLabel,filePrefix)

    % collect raw data, normalize by plasm if needed
    switch rawOrPlasmNorm
        case 'raw'
            yData = ec.clustT.(varName);    
        case 'plasmNorm'
            plasmVarName = strrep(varName,'clust_','plasm_');
            yData = ec.clustT.(varName)./ec.clustT.(plasmVarName);
            varName = [varName,'_plasmNorm'];
    end
    
    % setup labels
    if contains(varName,['C',num2str(channelIdx)])
        varNameForTitle = strrep(varName,...
            ['C',num2str(channelIdx)],...
            [' ',channelNamesForDisplays{channelIdx},' ']);
    else
        varNameForTitle = varName;
    end
    varNameForTitle = strrep(varNameForTitle,'_',' ');
    
    figTitle = ['hist of ',filePrefix,' ',varNameForTitle,' by Clusters'];
    
    plotArbitraryClusterVarHistogram(ec,yData,varName,figTitle,ecStagesToInclude,...
        idxData,conditionsOrder,condColors,...
        minVolume,maxVolume,histBinSize,...
        clustFolder,logTransform,histXLabel,filePrefix);    
end

%% plot a histogram of yData, arbitrary variable matching the cluster table size.
% plots a histogram for each condition, all in the same plot.
% also saves the figures and tables in the clustFolder

% INPUT
% ec: eggChamber object
% yData: variable - should be nCluster x 1 vector.
% varName: name of the variable to use in plot name and saved file names
% figTitle: title of the figure
% ecStagesToInclude: array of stages, e.g. [7,8] or 0:10 to include all
    % eggchambers if the staging wasnt performed.
% idxData: nCluster x 1  logical of which rows to include in the
    % computation
% conditionsOrder: order of plotting the conditions - e.g. [3,4,2,1] will
    % plot condition 3, then 4, then 2 then 1.
% condColors: colors assigned to each condition, Conditions x 3 array
% minVolume/maxVolume: min and maxVolume allowed for the cluster to be included
% hitBinSize: bin size of the histogram
% clustFolder: subfolder where tables/figs will be saved
% logTransform: whether to perform the histogram on a log10 of the data
% histXLabel: label for the x axis
function plotArbitraryClusterVarHistogram(ec,yData,varName,figTitle,ecStagesToInclude,...
    idxData,conditionsOrder,condColors,minVolume,maxVolume,histBinSize,...
    clustFolder,logTransform,histXLabel,filePrefix)
    
    logBase = 10;

    idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
            & idxData...
            & ec.clustT.clust_Volume >= minVolume...
            & ec.clustT.clust_Volume <= maxVolume;
    
    % apply log10 transformation if needed
    if logTransform
        [yData,~,~,warningMsg] = logTransformData(yData,idx0,logBase,varName);
    else 
        warningMsg = findWrongDataPoints(yData,idx0,'');
    end
    
    fh = figure('Name',figTitle);
    hold;
    for j=1:numel(ec.condIndices)
        curY = yData( idx0 ...
                        & ec.clustT.cond_Idx == ec.condIndices(conditionsOrder(j)) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY, min(curY(~isinf(curY))):histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        if numel(x) == 0
            disp(['cannot plot ',figTitle,' ',ec.conditionNames{conditionsOrder(j)},', check histogram bin size:']);
            disp(['   min Value: ',num2str(min(curY(~isinf(curY)))),'; max value: ',num2str(max(curY(~isinf(curY)))+histBinSize),'; histBinSize: ',num2str(histBinSize)]);
        else
            plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{conditionsOrder(j)});
            
            t = table(x',n'/sum(n),'VariableNames',{varName,'Count'});
            writetable(t,fullfile(clustFolder,['Hist',filePrefix,varName,'.txt']),'Delimiter','\t');
        end
        
    end
    if ~isempty(warningMsg{1})
        dim = [.5 .1 .3 .3];
        annotation('textbox',dim,'String',warningMsg,'FitBoxToText','on',...
            'Units','normalized','BackgroundColor','white');
    end
    legend show
    grid on
    xlabel(histXLabel);
    ylabel('Count');
    if logTransform
        xticklabels(10.^xticks); % label the ticks on the x-axis with the actual values, not the log-transformed ones
    end
    saveas(fh,fullfile(clustFolder,['Hist',filePrefix,varName,'.fig']));
    saveas(fh,fullfile(clustFolder,['Hist',filePrefix,varName,'.eps']),'epsc');
end

%%
function plotNucVarHistogramGvsS(ec,varName,ecStagesToInclude,...
    idxData,conditionsOrder,condColors,histBinSize,...
    channelNamesForDisplays,channelIdx,rawOrPlasmNorm,nucFolder,...
    logTransform,histXLabel,filePrefix,GvsSChannel,GvsSThresh,largeClusterMetricName)

    % collect raw data, normalize by plasm if needed
    switch rawOrPlasmNorm
        case 'raw'
            yData = ec.nucFullT.(varName); 
            yData(1:10)
        case 'plasmNorm'
            [~,~,baseName,suffix] = ec.splitVarName(varName);
            if contains(varName,['C',num2str(channelIdx)])
                varType = 'channel';
            else
                varType = 'geom';
            end
            plasmVarName = ec.buildVarName('plasm',channelIdx,baseName,suffix,varType);
            yData = ec.nucFullT.(varName)./ec.nucFullT.(plasmVarName);
            yData(1:10)
    end
    
    % setup labels
    if contains(varName,['C',num2str(channelIdx)])
        varNameForTitle = strrep(varName,...
            ['C',num2str(channelIdx)],...
            [' ',channelNamesForDisplays{channelIdx},' ']);
    else
        varNameForTitle = varName;
    end
    varNameForTitle = strrep(varNameForTitle,'_',' ');
    
    figTitle = ['hist of ',filePrefix,' ',varNameForTitle,' by Clusters'];

    plotArbitraryNucVarHistogramGvsS(ec,yData,varName,figTitle,ecStagesToInclude,...
        idxData,conditionsOrder,condColors,histBinSize,...
        nucFolder,logTransform,histXLabel,GvsSChannel,GvsSThresh,largeClusterMetricName);
end

%% log Transform data and generate a Warning message if some points are negative
function [logData,nWrong,nNeg,warningMsg,newVarName] = logTransformData(inData,idxData,logBase,varName)
    
    [warningMsg,idxWrong] = findWrongDataPoints(inData,idxData,[]);
    nWrong = sum(idxWrong);
    nTot = sum(idxData & (~idxWrong));
    nNeg = nTot - sum(idxData & (~idxWrong) & inData > 0);
    if nNeg>0
        curMsg = {['Warning: ',num2str(nNeg),'/',num2str(nTot),' data points'],...
                ' had negative values and couldn''t be plotted.'};
        if isempty(warningMsg{1})
            warningMsg = curMsg;
        else
            warningMsg = [warningMsg,...
                curMsg];
        end
    end

    logData = log(inData)/log(logBase);
    logData(inData<0) = NaN;
    
    %new variable name
    newVarName = [varName,'_Log',num2str(logBase)];
end

%% find number of NaNs and Infs in input Data and generate warning message
function [warningMsg,idxWrong] = findWrongDataPoints(inData,idxData,warningMsg)
    if isempty(warningMsg)
        warningMsg = {''};
    end
    idxWrong = idxData & (isnan(inData) | isinf(inData));
    nWrong = sum(idxWrong);
    if nWrong > 0
        if ~isempty(warningMsg{1})
            warningMsg = [warningMsg,...
                [['Warning: ',num2str(nWrong),'/',num2str(sum(idxData)),' data points'],...
                    ' had NaN/Inf values and couldn''t be plotted.']];
        else
            warningMsg = {['Warning: ',num2str(nWrong),'/',num2str(sum(idxData)),' data points'],...
                    ' had NaN/Inf values and couldn''t be plotted.'};
        end
    end
end

%%
function plotArbitraryNucVarHistogramGvsS(ec,yData,varName,figTitle,ecStagesToInclude,...
    idxData,conditionsOrder,condColors,histBinSize,...
    nucFolder,logTransform,histXLabel,GvsSChannel,GvsSThresh,largeClusterMetricName)
    
    logBase = 10;
    % plot distribution of the number of cluster per nucleus in G vs S
    
    gDist = cell(numel(ec.condIndices),1);
    sDist = cell(numel(ec.condIndices),1);
    fh = figure('Name',figTitle);
    hold on;
    
    % gray-ish colors for G phase
    grayColors = (4*0.7*ones(size(condColors))+ condColors)/5;
    
    idx0 = ismember(ec.nucFullT.eggChamber_Stage , ecStagesToInclude) ...
            & idxData;
    
    GvsSMetric = ec.nucFullT.(...
        ['nucAvgClust',largeClusterMetricName,'_clust_C',num2str(GvsSChannel),'Median_nucleoliSubtracted']) ...
        ./ ec.nucFullT.(['plasm_C',num2str(GvsSChannel),'Median_nucleoliSubtracted']);
    
    GvsSMetric = log10(GvsSMetric);
    
    % apply log10 transformation if needed
    if logTransform
        [yData,~,~,warningMsg] = logTransformData(yData,idx0,logBase,varName);
    else 
        warningMsg = findWrongDataPoints(yData,idx0,'');
    end

    % loop through conditions   
    xMin = min(yData(idx0 & ~isnan(yData) & ~isinf(yData)));
    xMax = max(yData(idx0 & ~isnan(yData) & ~isinf(yData)));
    for j=1:numel(ec.condIndices)
    
        curIdx = idx0 & ec.nucFullT.cond_Idx == conditionsOrder(j);
        % collect numbers of nuclei from G nuclei (value of metric < threshold or NaN)
        gDist{conditionsOrder(j)} = yData(curIdx ...
            & (isnan(GvsSMetric) | GvsSMetric < GvsSThresh(conditionsOrder(j)) ) );
        nTot = sum(curIdx);
    
        % generate & plot histogram
        [n,x] = hist(gDist{conditionsOrder(j)},xMin:histBinSize:xMax);
        nG = sum(n);
        percentG = nG/nTot*100;
        plot(x,n/sum(n),'Color',grayColors(j,:),'LineWidth',2,...
                    'DisplayName',[ec.conditionNames{conditionsOrder(j)},', G phase nuclei (',...
                    num2str(percentG,2),'% nuclei); mean = ',num2str(mean(gDist{conditionsOrder(j)},'omitnan'))]);
        % add data to table to be saved
        if j==1
            curT = table(x',n'/sum(n),'VariableNames',{'nHLBs',['Count_',ec.conditionNames{conditionsOrder(j)},'_G']});
        else
            curT = addvars(curT,n'/sum(n),'NewVariableNames',{['Count_',ec.conditionNames{conditionsOrder(j)},'_G']});
        end
    
        % collect numbers of nuclei from S nuclei (value of metric < threshold)
        sDist{conditionsOrder(j)} = yData(curIdx ...
            & GvsSMetric >= GvsSThresh(conditionsOrder(j))  );
         
        % generate & plot histogram
        [n,x] = hist(sDist{conditionsOrder(j)},xMin:histBinSize:xMax);
        nS =sum(n);
        percentS = nS/nTot*100;
        plot(x,n/sum(n),'Color',condColors(j,:),'LineWidth',2,...
                    'DisplayName',[ec.conditionNames{conditionsOrder(j)},', S phase nuclei (',...
                    num2str(percentS,2),'% nuclei); mean = ',num2str(mean(sDist{conditionsOrder(j)},'omitnan'))]);
        % add data to table to be saved
        curT = addvars(curT,n'/sum(n),'NewVariableNames',{['Count_',ec.conditionNames{conditionsOrder(j)},'_S']});
    end
    
    xlabel(histXLabel);
    ylabel('Count');
    grid on;
    legend show;

    if ~isempty(warningMsg{1})
        dim = [.5 .1 .3 .3];
        annotation('textbox',dim,'String',warningMsg,'FitBoxToText','on',...
            'Units','normalized','BackgroundColor','white');
    end
    
    % save figure and data
    figTitle = strrep(figTitle,' ','_');
    saveas(fh,fullfile(nucFolder,[figTitle,'.fig']));
    try
        saveas(fh,fullfile(nucFolder,[figTitle,'.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',figTitle,' - likely permission issue.']);
    end
    writetable(curT,fullfile(nucFolder,[figTitle,'.txt']),'Delimiter','\t');
end

%%
function scatterClusterVarChannelVsChannel(ec,varPrefix,varBaseName,varSuffix,...
    rawOrPlasmNorm,idxData,chX,chY,minVolume,maxVolume,ecStagesToInclude,conditionsOrder,...
    condColors,channelNamesForDisplays,alphaVal,qcFolder,filePrefix)

    % indices of points to plot
    idx0 = idxData...
        & ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
        & ec.clustT.eggChamber_Idx > 0 ...
        & ec.clustT.clust_Volume >= minVolume ...
        & ec.clustT.clust_Volume <= maxVolume;
    
    idx = cell(numel(ec.condIndices),1);
    for i=1:numel(ec.condIndices)
        idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(conditionsOrder(i));
    end
    
    fh = figure('Name',[filePrefix,' cluster Intensity, ',channelNamesForDisplays{chX}, ...
        ' (x) vs ',channelNamesForDisplays{chY},' (y) ',rawOrPlasmNorm]); 

    hold;
    varX = ec.buildVarName(varPrefix,chX,varBaseName,varSuffix,'channel');
    varY = ec.buildVarName(varPrefix,chY,varBaseName,varSuffix,'channel');

    for i=1:numel(ec.condIndices)
        switch rawOrPlasmNorm
            case 'raw'
                xData = ec.clustT.(varX)(idx{i});
                yData = ec.clustT.(varY)(idx{i});
            case 'plasmNorm'
                varDenomX = ec.buildVarName('plasm',chX,varBaseName,varSuffix,'channel');
                varDenomY = ec.buildVarName('plasm',chY,varBaseName,varSuffix,'channel');
                xData = ec.clustT.(varX)(idx{i})./ec.clustT.(varDenomX)(idx{i});
                yData = ec.clustT.(varY)(idx{i})./ec.clustT.(varDenomY)(idx{i});
        end

        scatter(xData,yData,...
            'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',...
            condColors(i,:),'DisplayName',ec.conditionNames{conditionsOrder(i)});
    end 
    
    alpha(alphaVal);
    xlabel(strrep([varX,' ',rawOrPlasmNorm],'_',' '));
    ylabel(strrep([varY,' ',rawOrPlasmNorm],'_',' '));
    grid on;
    legend show;

    saveas(fh,fullfile(qcFolder,...
        [filePrefix,'clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'.fig']));

    try
        saveas(fh,fullfile(qcFolder,...
        [filePrefix,'clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',...
            ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'.eps'],...
            ' - likely permission issue.']);
    end
end

