fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img';

%% set path for subfunctions

addpath('subfunctions/');

%%

% intitialize eggChamberDataFolder objec - this collecs the locations of
% the files associated with the experiment and stores them by condition and sample
ec = eggChamberDataSet(fijiOutFolder);

%% load nucleus-wide raw data
% load raw data from nuclei segmentation - builds a table ec.nucFullT that combines
% nuclei from all samples (one row per nucleus) and compiles dozens of
% metrics of size and intensity in every channel. (creates nucT as a backup copy)
ec.loadAllEggChamberNucleiData;
%%
% cleans up the ec.nucT table by:
    % removing useless metrics
    % removing data from the dummy channel used as a marker of eggChamberID
    % sorting the columns so that the sample/eggchamber/nucleus info is in
    % the first columns
    % ordering the rows by eggChamberID.
% ec.nucFullT retains all metrics.
%ec.streamLineTable('nuc');

%% load cluster data
% populates the table ec.clustT
% also adds nucleoli and nucleoplasm data to ec.nucT and ec.nucFullT tables
% does NOT add summary cluster metrics to nucT and nucFullT
ec.loadAllClusterData();

%% background correct nuclei intensity
ec.backgroundCorrectNucIntensity;

%% add nuc stats to cluster table
% takes metrics from nucT and copies them as extra variables of the cluster table
ec.clustT = ec.addNucStatsToClustTable();

%% ackground correct clusters intensity
ec.backgroundCorrectClustIntensity();

%% add cluster data to nuc table
% takes metrics from clustT, averages them over each nucleus
% and copies them as extra variables of the nuc table
ec.addAverageClusterStatsToNucTable();

%% QC block
figure; hold;
[~,idx] = sort(ec.nucFullT.nucSphericity,'descend');
plot(zscore(ec.nucFullT.nucSphericity(idx)),'DisplayName','sphericity');
plot(zscore(ec.nucFullT.nucMeanBreadth(idx)),'DisplayName','Breadth');
plot(zscore(ec.nucFullT.nucSurfaceArea(idx)),'DisplayName','surf Area');
plot(zscore(tJoin.numClusters(idx)),'DisplayName','num Clusters');
plot(zscore(ec.nucFullT.nucVolume(idx)),'DisplayName','Volume');

%%

figure; hold;
idx = ismember(ec.nucFullT.eggChamberStage,6);
plot(ec.nucFullT.nucSphericity(idx),ec.nucFullT.nucVolume(idx),'o','MarkerFaceColor','auto');
idx = ismember(ec.nucFullT.eggChamberStage,7);
plot(ec.nucFullT.nucSphericity(idx),ec.nucFullT.nucVolume(idx),'o','MarkerFaceColor','auto');
idx = ismember(ec.nucFullT.eggChamberStage,8);
plot(ec.nucFullT.nucSphericity(idx),ec.nucFullT.nucVolume(idx),'o','MarkerFaceColor','auto');

figure; hold;
idx = ismember(ec.nucFullT.eggChamberStage,6);
plot(modzscore(ec.nucFullT.nucSphericity(idx)),modzscore(ec.nucFullT.nucVolume(idx)),'o','MarkerFaceColor','auto');
idx = ismember(ec.nucFullT.eggChamberStage,7);
plot(modzscore(ec.nucFullT.nucSphericity(idx)),modzscore(ec.nucFullT.nucVolume(idx)),'o','MarkerFaceColor','auto');
idx = ismember(ec.nucFullT.eggChamberStage,8);
plot(modzscore(ec.nucFullT.nucSphericity(idx)),modzscore(ec.nucFullT.nucVolume(idx)),'o','MarkerFaceColor','auto');

%% generate summary statistics per egg chamber
sumT = ec.generateEggChamberSummaryTable();

%% compute background-correced nuclei intensity values (eggChamber and whole image)
ec.backgroundCorrectNucIntensity();

%% scatter plot by sample
% plot where all nuclei from each sample are grouped separately.
% not terribly useful unless the egg chambers havent been segmented yet.
ec.scatterPlotNucleiMetricBySample('nuc',1,'Volume','');
ec.scatterPlotNucleiMetricBySample('nuc',1,'Mean','wholeImgCorr');

%% scatter plot by egg chamber
% plot where all nuclei from an egg chamber are grouped separately.
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Volume','',[6,7,8]); % nuclear Volume

% nuclear mean intensity in DAPI channel (raw)
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Mean','raw',[6,7,8]); 

% nuclear mean intensity in DAPI channel 
% (correced for background estimated across the egg chambers in the image)
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Mean','wholeImgCorr','all'); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Mean','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Mean','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Mean','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Mean','wholeImgCorr',[6,7,8]); 

ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','wholeImgCorr','all'); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Median','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Median','wholeImgCorr',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Median','wholeImgCorr',[6,7,8]); 

%% plot histograms

