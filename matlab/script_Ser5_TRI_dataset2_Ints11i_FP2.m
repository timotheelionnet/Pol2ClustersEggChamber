% enter here the path of the folder where the output of the Fiji scripts is saved
fijiOutFolder = '/Volumes/Lionnet1/20230812-matlab/20230812-IntS11ivsCtrl_Ser5ph-out';
fijiOutFolder = '/Volumes/lionnt01lab/lionnt01labspace/Feiyue_Tim/20230812-Ser5ph dataset with IntS11i vs control, FP vs control/20230812-IntS11ivsCtrl_Ser5ph-out';

% entere here the path where to re-load the data from a previous Matlab run
% (only used if useMatlabRatherThanFiji = 1;)
matlabInFolder = '/Volumes/lionnt01lab/lionnt01labspace/Feiyue_Tim/20230812-Ser5ph dataset with IntS11i vs control, FP vs control/20230812-IntS11ivsCtrl_Ser5ph-matlabout';

% entere here the path where to save the data 
matlabOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/Ints11i_Ser5ph_TRI/matlabOut2';

% if the data has already been uploaded in Matlab and saved in the
% matlabOutFolder, set this Flag to 1 (faster)
useMatlabRatherThanFiji = 0;

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

%% scatter plot by egg chamber: Nucleus Volume (qc)
% plot where all nuclei from an egg chamber are grouped separately.

% all stages
[fullTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
    'nuc',1,'Volume','','nucVolume',nucleiInGoodEggChambers,{},...
    conditionsOrder,ecStagesToInclude,1,'useMean',0.1);

saveNucDataFromPlot(fh,fullTable,avgEcTable,avgCondTable, qcFolder,'nucVolume');

%% scatter plot nuclei intensity by egg chamber (qc): 

% raw:
for i=1:numel(channelNamesForDisplays)
    [nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
        'nuc',i,'Mean','raw',...
        ['meanNucInt',channelNamesForDisplays{i},'_raw'],true(size(ec.nucFullT,1),1),...
        {},conditionsOrder,ecStagesToInclude,1,'useMean',0.3);

    saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, qcFolder,...
        ['meanNucInt',channelNamesForDisplays{i},'_raw']);
end

% eggChamber-corrected
for i=1:numel(channelNamesForDisplays)

    [nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
        'nuc',i,'Mean','eggChamberSubtracted',...
        ['meanNucInt',channelNamesForDisplays{i},'_eggChamberSubtracted'],true(size(ec.nucFullT,1),1),...
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
        ['meanPlasmInt',channelNamesForDisplays{i},'_nucleoliSub'],true(size(ec.nucFullT,1),1),...
        {},conditionsOrder,ecStagesToInclude,1,'useMean',0.3);

    saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, nucFolder,...
        ['meanPlasmInt',channelNamesForDisplays{i},'_nucleoliSub']);
end

%% cluster intensity as a function of cluster volume, raw (qc)
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0;

for i=1:numel(channelNamesForDisplays)
    fh = figure('Name',...
        ['Cluster Volume vs ',channelNamesForDisplays{i},' intensity']);

    s = scatter(ec.clustT.clust_Volume(idx0),...
        ec.clustT.(['clust_C',num2str(i),'Mean_plasmCorr'])(idx0),10,...
        ec.clustT.cond_Idx(idx0),'o','filled');
    alpha(s,0.1);
    colormap(condColors);
    colorbar('Ticks',1:numel(ec.condIndices),'Ticklabels',ec.conditionNames);
    xlabel('Volume (um^3)');
    ylabel(['clust C',num2str(i),'Mean plasmCorr']);
    set(gca,'xscale','log');
    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'.fig']));
    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'.eps']),'epsc');
end

%% cluster intensity as a function of cluster volume, normalized by nucleoplasm intensity
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0;

for i=1:numel(channelNamesForDisplays)
    fh = figure('Name',...
        ['Cluster Volume vs ',channelNamesForDisplays{i},' intensity, nucleoplasm-normalized']);

    s = scatter(ec.clustT.clust_Volume(idx0),...
        log(ec.clustT.(['clust_C',num2str(i),'Mean_plasmCorr'])(idx0)...
        ./ec.clustT.(['plasm_C',num2str(i),'Mean_eggChamberSubtracted'])(idx0))/log(2),...
        10,ec.clustT.cond_Idx(idx0),'o','filled');

    alpha(s,0.1);
    colormap(condColors);
    colorbar('Ticks',1:numel(ec.condIndices),'Ticklabels',ec.conditionNames);
    xlabel('Volume (um^3)');
    ylabel(['log2 clust C',num2str(i),'Mean plasmCorr']);
    set(gca,'xscale','log');
    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'_plasmNorm.fig']));
    saveas(fh,fullfile(qcFolder,['clustVolvs',channelNamesForDisplays{i},'_plasmNorm.eps']),'epsc');
end

%% scatter plot by egg chamber: Number of large clusters per nucleus, Cluster volume
minVolume = hlbMinVol;
maxVolume = Inf;

% plot number of large clusters per nucleus
[nucTable,avgEcTable,avgCondTable,fh] = ec.scatterPlotNucTableMetricByEggChamber( ...
    'nuc',1,'NumClustersMinVol','','nLargeClusters',true(size(ec.nucT,1),1),...
    {},conditionsOrder,ecStagesToInclude,1,'useMean',0.2);

saveNucDataFromPlot(fh,nucTable,avgEcTable,avgCondTable, nucFolder,'nClusters');

% plot cluster volume for large clusters
[clustTable,avgNucClustTable,avgEcTable,avgCondTable,fhFull] = ec.scatterPlotClustTableMetricByEggChamber( ...
    'clust',1,'Volume','','clusterVolume',true(size(ec.clustT,1),1),...
    {},conditionsOrder,ecStagesToInclude,minVolume,maxVolume,1,'useMean',0.1);

% same, but nucleus average
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'mean_clusterVolume',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'mean clusterVolume by Nuclei',0.1);

% save tables
saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcTable,avgCondTable, clustFolder,'clustVolume');

%% scatter plot by egg chamber: Number of small clusters per nucleus, Cluster volume
minVolume = hlbMinVol/20;
maxVolume = hlbMinVol/5;

% plot cluster volume for small clusters
[clustTable,avgNucClustTable,avgEcTable,avgCondTable,fhFull] = ec.scatterPlotClustTableMetricByEggChamber( ...
    'clust',1,'Volume','','clusterVolume',true(size(ec.clustT,1),1),...
    {},conditionsOrder,ecStagesToInclude,minVolume,maxVolume,1,'useMean',0.01);

% same, but nucleus average
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'mean_clusterVolume',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'mean clusterVolume by Nuclei',0.05);

% plot the number of small clusters per nucleus
fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,'nClusters',...
    true(size(avgNucClustTable,1),1),...
    conditionsOrder,ecStagesToInclude,'number of Small Clusters per Nucleus',0.05);

% save tables
saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcTable,avgCondTable, clustFolder,'smallClustVolume');

%% cluster intensity, HLBs
minVolume =hlbMinVol;
maxVolume = Inf;

for i=1:numel(channelNamesForDisplays)
    curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustTableMetricByEggChamber(...
        'clust',i,'Median','nucleoliSubtracted',curVarName,...
        true(size(ec.clustT,1),1),{},conditionsOrder,ecStagesToInclude,...
        minVolume,maxVolume,1,'useMean',0.05);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['clustInt',channelNamesForDisplays{i},'_nucleoliSub']);
end

%% histogram of HLB volumes by condition 
minVolume =hlbMinVol;
maxVolume = Inf;
histBinSize = 0.03;

% collect all intensities of the current channel
yData = log(ec.clustT.('clust_Volume'))/log(10);

fh = figure('Name','histogram of large cluster Volume condition.');
hold;
for j=1:numel(ec.condIndices)
    curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                    & ec.clustT.clust_Volume >= minVolume...
                    & ec.clustT.clust_Volume <= maxVolume);
    [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
    plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

    t = table(x',n'/sum(n),'VariableNames',{'Clust_Volume_Log10','Count'});
    writetable(t,fullfile(clustFolder,'HistClustVolumeLog10.txt'),'Delimiter','\t');
    
end
legend show
grid on
xlabel('Cluster Volume (um^3, Log scale)');
ylabel('Count');
xticklabels(10.^xticks); % label the ticks on the x-axis with the actual values, not the log-transformed ones

saveas(fh,fullfile(clustFolder,'HistClustVolume.fig'));
saveas(fh,fullfile(clustFolder,'HistClustVolume.eps'),'epsc');
    

%% cluster intensity for SMALL clusters (nucleoli subtracted)
minVolume =0;
maxVolume = hlbMinVol/5;

for i=1:numel(channelNamesForDisplays)
   curVarName = ['SmallClustInt_',channelNamesForDisplays{i},'_nucleoliSubtracted'];
    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustTableMetricByEggChamber(...
        'clust',i,'Median','nucleoliSubtracted',curVarName,...
        true(size(ec.clustT,1),1),{},conditionsOrder,ecStagesToInclude,...
        minVolume,maxVolume,1,'useMean',0.01);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['smallClustInt',channelNamesForDisplays{i},'_nucleoliSub']);
end

%% cluster intensity, HLBs (nucleoli subtracted, relative to plasm levels)
minVolume =hlbMinVol;
maxVolume = Inf;

for i=1:numel(channelNamesForDisplays)
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);
    
    curVarName = ['clustInt_',channelNamesForDisplays{i},'_nucleoliSubtr_plasmNorm'];

    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustArbitraryMetricByEggChamber(yData,curVarName,...
        true(size(ec.clustT,1),1),{},conditionsOrder,ecStagesToInclude,...
        minVolume,maxVolume,1,'useMean',0.2);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['clustInt_',channelNamesForDisplays{i},'_plasmNorm']);
end


%% plot number of clusters vs cluster signal



%% cluster intensity for SMALL clusters (nucleoli subtracted, relative to plasm levels)
minVolume = 0;
maxVolume = hlbMinVol/5;

for i=1:numel(channelNamesForDisplays)
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);

    curVarName = ['smallClustInt_',channelNamesForDisplays{i},'_nucleoliSubtr_plasmNorm'];

    [clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable,fhFull] = ...
        ec.scatterPlotClustArbitraryMetricByEggChamber(...
        yData,curVarName,...
        true(size(ec.clustT,1),1),{},conditionsOrder,ecStagesToInclude,...
        minVolume,maxVolume,1,'useMean',0.01);

    fhNuc = ec.plotClusterMetricByNucleus(avgNucClustTable,...
        ['mean_',curVarName],true(size(avgNucClustTable,1),1),...
        conditionsOrder,ecStagesToInclude,['mean ',curVarName,' by Nuclei'],0.2);

    saveClustDataFromPlot(fhFull,fhNuc,clustTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['smallClustInt_',channelNamesForDisplays{i},'_plasmNorm']);
end

%% histogram of HLB cluster intensities (relative to plasm levels) by condition
minVolume =hlbMinVol;
maxVolume = Inf;
histBinSize = 0.2;

for i=1:numel(channelNamesForDisplays)
    % collect all intensities of the current channel
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);

    fh = figure('Name',['histogram of cluster ',channelNamesForDisplays{i},' Int by condition.']);
    hold;
    for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n'/sum(n),'VariableNames',{['Clust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm'],'Count'});
        writetable(t,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.txt']),'Delimiter','\t');
        
    end
    legend show
    grid on
    xlabel(['Cluster ',channelNamesForDisplays{i},' intensity, plasm-norm']);
    ylabel('Count');
    saveas(fh,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.fig']));
    saveas(fh,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.eps']),'epsc');
    
end

%% histogram of SMALL cluster intensities (relative to plasm levels) by condition
minVolume =0;
maxVolume = hlbMinVol/5;
histBinSize = 0.2;

for i=1:numel(channelNamesForDisplays)
    % collect all intensities of the current channel
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);

    fh = figure('Name',['histogram of SMALL cluster ',channelNamesForDisplays{i},' Int by condition.']);
    hold;
    for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n'/sum(n),'VariableNames',{['SmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm'],'Count'});
        writetable(t,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.txt']),'Delimiter','\t');
    end
    legend show
    grid on
    xlabel(['Small Cluster ',channelNamesForDisplays{i},' intensity, plasm-norm']);
    ylabel('Count');
    saveas(fh,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.fig']));
    saveas(fh,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntPlasmNorm.eps']),'epsc');
    
end

%% histogram of HLB cluster intensities (raw) by condition
minVolume =hlbMinVol;
maxVolume = Inf;
histBinSize = 100;

for i=1:numel(channelNamesForDisplays)
    % collect all intensities of the current channel
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']);

    fh = figure('Name',['histogram of raw HLB ',channelNamesForDisplays{i},' Int by condition (bg-corr).']);
    hold;
    for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n'/sum(n),'VariableNames',{['Clust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr'],'Count'});
        writetable(t,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.txt']),'Delimiter','\t');
        
    end
    legend show
    xlabel(['Cluster ',channelNamesForDisplays{i},' intensity, nucleoli-subtracted']);
    ylabel('Count');
    saveas(fh,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.fig']));
    saveas(fh,fullfile(clustFolder,['HistClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.eps']),'epsc');
    
end

%% histogram of SMALL cluster intensities (raw) by condition
minVolume =0;
maxVolume = hlbMinVol/5;
histBinSize = 50;

for i=1:numel(channelNamesForDisplays)
    % collect all intensities of the current channel
    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']);

    fh = figure('Name',['histogram of raw SMALL cluster ',channelNamesForDisplays{i},' Int by condition (bg corr).']);
    hold;
    for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n/sum(n),'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n'/sum(n),'VariableNames',{['SmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr'],'Count'});
        writetable(t,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.txt']),'Delimiter','\t');
    end
    legend show
    xlabel(['Small Cluster ',channelNamesForDisplays{i},' intensity, nucleoli-subtracted']);
    ylabel('Count');
    saveas(fh,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.fig']));
    saveas(fh,fullfile(clustFolder,['HistSmallClust',channelNamesForDisplays{i},'_',ec.conditionNames{j},'IntBgCorr.eps']),'epsc');
    
end

%% Cluster intensity, channel vs channel (nucleoli subtracted)
chX = 3;  % channel for X axis
chY = 4;  % channel for Y axis
minVolume = hlbMinVol;
maxVolume = Inf;

% indices of points to plot
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0 ...
    & ec.clustT.clust_Volume >= minVolume ...
    & ec.clustT.clust_Volume <= maxVolume;

idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

fh = figure('Name',['cluster Intensity, ',channelNamesForDisplays{chX},' (x) vs ',channelNamesForDisplays{chY},' (y)']); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i}),...
        ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 

alpha(0.3);
xlabel(['C',num2str(chX),'Median nucleoliSubtracted']);
ylabel(['C',num2str(chY),'Median nucleoliSubtracted']);
grid on;
legend show;

saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'.fig']));
saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'.eps']),'epsc');

%% Cluster intensity, channel vs channel (relative to nuclei levels)
chX = 4; % channel for X axis
chY = 3; % channel for Y axis
minVolume = hlbMinVol;
maxVolume = Inf;

% indices of points to plot
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0 ...
    & ec.clustT.clust_Volume >= minVolume ...
    & ec.clustT.clust_Volume <= maxVolume;

idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

fh = figure('Name',['cluster Intensity, ',channelNamesForDisplays{chX},' (x) vs ',channelNamesForDisplays{chY},' (y), nucleoplasm-normalized']); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i}),...
        ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 

alpha(0.3);
xlabel(['C',num2str(chX),'Median nucleoliSubtracted, plasmNorm']);
ylabel(['C',num2str(chY),'Median nucleoliSubtracted, plasmNorm']);
grid on;
legend show;

saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalize.fig']));
saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalized.eps']),'epsc');

clear idx idx0

%% Cluster intensity, channel vs channel (relative to nuclei levels)
chX = 2; % channel for X axis
chY = 3; % channel for Y axis
minVolume = hlbMinVol;
maxVolume = Inf;

% indices of points to plot
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0 ...
    & ec.clustT.clust_Volume >= minVolume ...
    & ec.clustT.clust_Volume <= maxVolume;

idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

fh = figure('Name',['cluster Intensity, ',channelNamesForDisplays{chX},' (x) vs ',channelNamesForDisplays{chY},' (y), nucleoplasm-normalized']); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i}),...
        ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 

alpha(0.3);
xlabel(['C',num2str(chX),'Median nucleoliSubtracted, plasmNorm']);
ylabel(['C',num2str(chY),'Median nucleoliSubtracted, plasmNorm']);
grid on;
legend show;

saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalize.fig']));
saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalized.eps']),'epsc');

clear idx idx0

%% Small Cluster intensity, channel vs channel (relative to nuclei levels)
chX = 4; % channel for X axis
chY = 3; % channel for Y axis
minVolume = 0;
maxVolume = hlbMinVol/5;

% indices of points to plot
idx0 = ismember(ec.clustT.eggChamber_Stage , ecStagesToInclude) ...
    & ec.clustT.eggChamber_Idx > 0 ...
    & ec.clustT.clust_Volume >= minVolume ...
    & ec.clustT.clust_Volume <= maxVolume;

idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

fh = figure('Name',['Small cluster Intensity, ',channelNamesForDisplays{chX},' (x) vs ',channelNamesForDisplays{chY},' (y), nucleoplasm-normalized']); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted'])(idx{i}),...
        ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i})...
        ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted'])(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 

alpha(0.01);
xlabel(['C',num2str(chX),'Median nucleoliSubtracted, plasmNorm']);
ylabel(['C',num2str(chY),'Median nucleoliSubtracted, plasmNorm']);
grid on;
legend show;

saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalize.fig']));
saveas(fh,fullfile(qcFolder,...
    ['clustInt',channelNamesForDisplays{chX},'_vs_',channelNamesForDisplays{chY},'_plasmNormalized.eps']),'epsc');

clear idx idx0
%% HLB chY/chX ratio, by sample - Normalized to nuclear levels (nucleoli Subtracted, normalized to nucleoplasm levels)
chY = 3; % numerator channel
chX = 4; % denominator channel

minVolume = hlbMinVol; % minimum cluster Volume
maxVolume = Inf; % max cluster Volume

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
minVolume =hlbMinVol;
maxVolume = Inf;
histBinSize = 0.1;

chY = 3; % numerator channel
chX = 4; % denominator channel

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

fh = figure('Name',['Hist HLB ',channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX},'normalized by nucleoplasm']);
hold;
for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n,'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n','VariableNames',{['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'],'Count'});
        writetable(t,fullfile(clustFolder,['HistClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.txt']),'Delimiter','\t');
end

legend show
xlabel(['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm']);
ylabel('Count');
saveas(fh,fullfile(clustFolder,['HistClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.fig']));
saveas(fh,fullfile(clustFolder,['HistClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.eps']),'epsc');

%% histogram of SMALL cluster chY/chX ratio (relative to plasm levels) by condition
minVolume = 0 ;
maxVolume = hlbMinVol/5;
histBinSize = 0.1;

chY = 3; % numerator channel
chX = 4; % denominator channel

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

fh = figure('Name',['small cluster ',channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX},' normalized by nucleoplasm']);
hold;
for j=1:numel(ec.condIndices)
        curY = yData( ec.clustT.cond_Idx ==ec.condIndices(j) ...
                        & ec.clustT.clust_Volume >= minVolume...
                        & ec.clustT.clust_Volume <= maxVolume);
        [n,x] = hist(curY,0:histBinSize:(max(curY(~isinf(curY)))+histBinSize));
        plot(x,n,'-','Color',condColors(j,:),'LineWidth',2,'DisplayName',ec.conditionNames{j});

        t = table(x',n','VariableNames',{['smallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm'],'Count'});
        writetable(t,fullfile(clustFolder,['HistSmallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.txt']),'Delimiter','\t');
end

legend show
xlabel(['smallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm']);
ylabel('Count');
xlim([0 10]);
saveas(fh,fullfile(clustFolder,['HistSmallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.fig']));
saveas(fh,fullfile(clustFolder,['HistSmallClustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm.eps']),'epsc');

%% pl0tting MPM2+ and MPM2- HLB metrics side by side
chRef = 2; % channel used to startify the data in positive vs negative
yMax1 = 2.5; % max value to be called a negative
yMin2 = 5; % min value to be called a positive

stratificationName = 'MPM2positive'; %name of stratification variable
idx1StratificationVal = 0; % negative
idx2StratificationVal = 1; % positive

minVolume = hlbMinVol; % minimum cluster Volume
maxVolume = Inf; % max cluster Volume

yRef = ec.clustT.(['clust_C',num2str(chRef),'Median_nucleoliSubtracted']) ...
    ./ec.clustT.(['plasm_C',num2str(chRef),'Median_nucleoliSubtracted']);

idxData1 = yRef < yMax1;
idxData2 = yRef > yMin2;

% first category (idxData1) is gray, second (idxData2) is in the default colors
dualColorMap = {0.7*ones(size(condColors));condColors};

for i=1:numel(channelNamesForDisplays)

    yData = ec.clustT.(['clust_C',num2str(i),'Median_nucleoliSubtracted']) ...
        ./ec.clustT.(['plasm_C',num2str(i),'Median_nucleoliSubtracted']);

    [clustTable1,avgEcClustTable1,avgCondClustTable1,...
    clustTable2,avgEcClustTable2,avgCondClustTable2,fh] = ...
                ec.scatterPlotClustDualArbitraryMetricByEggChamber(...
                yData,['clustInt_',channelNamesForDisplays{i},'_plasmNorm'],idxData1,idxData2,...
                {},conditionsOrder,ecStagesToInclude,...
                minVolume,maxVolume,1,'useMean',0.3,dualColorMap);
    ylim('auto');

    clustTable1 = addvars(clustTable1,idx1StratificationVal*ones(size(clustTable1,1),1),'newVariableNames',{stratificationName});
    avgEcClustTable1 = addvars(avgEcClustTable1,idx1StratificationVal*ones(size(avgEcClustTable1,1),1),'newVariableNames',{stratificationName});
    avgCondClustTable1 = addvars(avgCondClustTable1,idx1StratificationVal*ones(size(avgCondClustTable1,1),1),'newVariableNames',{stratificationName});
    
    clustTable2 = addvars(clustTable2,idx2StratificationVal*ones(size(clustTable2,1),1),'newVariableNames',{stratificationName});
    avgEcClustTable2 = addvars(avgEcClustTable2,idx2StratificationVal*ones(size(avgEcClustTable2,1),1),'newVariableNames',{stratificationName});
    avgCondClustTable2 = addvars(avgCondClustTable2,idx2StratificationVal*ones(size(avgCondClustTable2,1),1),'newVariableNames',{stratificationName});
    
    clustTable = [clustTable1;clustTable2];
    avgEcClustTable = [avgEcClustTable1;avgEcClustTable2];
    avgCondClustTable = [avgCondClustTable1;avgCondClustTable2];
    
    saveNucDataFromPlot(fh,clustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
            ['clustInt_',channelNamesForDisplays{i},'_plasmNorm_MPM2strat']);

end

%% pltting MPM2+ and MPM2- Ser5ph/Pol II side by side
chRef = 2; % channel used to startify the data in positive vs negative
yMax1 = 2.5; % max value to be called a negative
yMin2 = 5; % min value to be called a positive
stratificationName = 'MPM2positive';
idx1StratificationVal = 0;
idx2StratificationVal = 1;

chY = 3;
chX = 4;

minVolume = hlbMinVol; % minimum cluster Volume
maxVolume = Inf; % max cluster Volume

yRef = ec.clustT.(['clust_C',num2str(chRef),'Median_nucleoliSubtracted']) ...
    ./ec.clustT.(['plasm_C',num2str(chRef),'Median_nucleoliSubtracted']);

idxData1 = yRef < yMax1;
idxData2 = yRef > yMin2;

% first category (idxData1) is gray, second (idxData2) is in the default colors
dualColorMap = {0.7*ones(size(condColors));condColors};

numData = ec.clustT.(['clust_C',num2str(chY),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chY),'Median_nucleoliSubtracted']);

denomData = ec.clustT.(['clust_C',num2str(chX),'Median_nucleoliSubtracted'])...
    ./ec.clustT.(['plasm_C',num2str(chX),'Median_nucleoliSubtracted']); 

yData = numData./denomData;

[clustTable1,avgEcClustTable1,avgCondClustTable1,...
    clustTable2,avgEcClustTable2,avgCondClustTable2,fh] = ...
    ec.scatterPlotClustDualArbitraryMetricByEggChamber(...
        yData,[channelNamesForDisplays{chY},'/',channelNamesForDisplays{chX}],idxData1,idxData2,...
        {},conditionsOrder,ecStagesToInclude,...
        minVolume,maxVolume,1,'useMean',0.3,dualColorMap);

clustTable1 = addvars(clustTable1,idx1StratificationVal*ones(size(clustTable1,1),1),'newVariableNames',{stratificationName});
avgEcClustTable1 = addvars(avgEcClustTable1,idx1StratificationVal*ones(size(avgEcClustTable1,1),1),'newVariableNames',{stratificationName});
avgCondClustTable1 = addvars(avgCondClustTable1,idx1StratificationVal*ones(size(avgCondClustTable1,1),1),'newVariableNames',{stratificationName});

clustTable2 = addvars(clustTable2,idx2StratificationVal*ones(size(clustTable2,1),1),'newVariableNames',{stratificationName});
avgEcClustTable2 = addvars(avgEcClustTable2,idx2StratificationVal*ones(size(avgEcClustTable2,1),1),'newVariableNames',{stratificationName});
avgCondClustTable2 = addvars(avgCondClustTable2,idx2StratificationVal*ones(size(avgCondClustTable2,1),1),'newVariableNames',{stratificationName});

clustTable = [clustTable1;clustTable2];
avgEcClustTable = [avgEcClustTable1;avgEcClustTable2];
avgCondClustTable = [avgCondClustTable1;avgCondClustTable2];

saveNucDataFromPlot(fh,clustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
        ['clustInt',channelNamesForDisplays{chX},'-',channelNamesForDisplays{chY},'ratio_plasmNorm_MPM2strat']);

%% pltting MPM2+ and MPM2- HLB Volume side by side
chRef = 2; % channel used to startify the data in positive vs negative
yMax1 = 1; % max value to be called a negative
yMin2 = 5; % min value to be called a positive

stratificationName = 'MPM2positive'; %name of stratification variable
idx1StratificationVal = 0; % negative
idx2StratificationVal = 1; % positive

minVolume = hlbMinVol; % minimum cluster Volume
maxVolume = Inf; % max cluster Volume

yRef = ec.clustT.(['clust_C',num2str(chRef),'Median_nucleoliSubtracted']) ...
    ./ec.clustT.(['plasm_C',num2str(chRef),'Median_nucleoliSubtracted']);

idxData1 = yRef < yMax1;
idxData2 = yRef > yMin2;

% first category (idxData1) is gray, second (idxData2) is in the default colors
dualColorMap = {0.7*ones(size(condColors));condColors};

yData = ec.clustT.clust_Volume;

[clustTable1,avgEcClustTable1,avgCondClustTable1,...
clustTable2,avgEcClustTable2,avgCondClustTable2,fh] = ...
            ec.scatterPlotClustDualArbitraryMetricByEggChamber(...
            yData,'clust_Volume',idxData1,idxData2,...
            {},conditionsOrder,ecStagesToInclude,...
            minVolume,maxVolume,1,'useMean',0.3,dualColorMap);
ylim('auto');

% clustTable1 = addvars(clustTable1,idx1StratificationVal*ones(size(clustTable1,1),1),'newVariableNames',{stratificationName});
% avgEcClustTable1 = addvars(avgEcClustTable1,idx1StratificationVal*ones(size(avgEcClustTable1,1),1),'newVariableNames',{stratificationName});
% avgCondClustTable1 = addvars(avgCondClustTable1,idx1StratificationVal*ones(size(avgCondClustTable1,1),1),'newVariableNames',{stratificationName});
% 
% clustTable2 = addvars(clustTable2,idx2StratificationVal*ones(size(clustTable2,1),1),'newVariableNames',{stratificationName});
% avgEcClustTable2 = addvars(avgEcClustTable2,idx2StratificationVal*ones(size(avgEcClustTable2,1),1),'newVariableNames',{stratificationName});
% avgCondClustTable2 = addvars(avgCondClustTable2,idx2StratificationVal*ones(size(avgCondClustTable2,1),1),'newVariableNames',{stratificationName});
% 
% clustTable = [clustTable1;clustTable2];
% avgEcClustTable = [avgEcClustTable1;avgEcClustTable2];
% avgCondClustTable = [avgCondClustTable1;avgCondClustTable2];
% 
% saveNucDataFromPlot(fh,clustTable,avgEcClustTable,avgCondClustTable, clustFolder,...
%         ['clustInt_',channelNamesForDisplays{i},'_plasmNorm_MPM2strat']);

