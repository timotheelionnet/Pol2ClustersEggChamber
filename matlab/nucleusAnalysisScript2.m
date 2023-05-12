fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img';

%% set path for subfunctions

addpath('subfunctions/');

%%

% intitialize eggChamberDataFolder object - this collects the locations of
% the files associated with the experiment and stores them by condition and sample
ec = eggChamberDataFolder(fijiOutFolder);

%% load nucleus-wide raw data
% load raw data from nuclei segmentation - builds a table t that combines
% nuclei from all samples (one row per nucleus) and compiles dozens of
% metrics of size and intensity in every channel
ect = ec.loadAllEggChamberData();

% cleans up the table by:
    % removing useless metrics
    % removing data from the dummy channel used as a marker of eggChamberID
    % sorting the columns so that the sample/eggchamber/nucleus info is in
    % the first columns
    % ordering the rows by eggChamberID.
ect.streamLineTable;

%% generate summary statistics per egg chamber
sumT = ect.generateEggChamberSummaryTable();

%% compute background-corrected nuclei intensity values (eggChamber and whole image)
ect.backgroundCorrectNucIntensity;

%% scatter plot by sample
% plot where all nuclei from each sample are grouped separately.
% not terribly useful unless the egg chambers havent been segmented yet.
ect.scatterPlotMetricBySample('nuc',1,'Volume','');
ect.scatterPlotMetricBySample('nuc',1,'Mean','wholeImgCorr');

%% scatter plot by egg chamber
% plot where all nuclei from an egg chamber are grouped separately.
ect.scatterPlotMetricByEggChamber('nuc',1,'Volume','',[6,7,8]); % nuclear Volume

% nuclear mean intensity in DAPI channel (raw)
ect.scatterPlotMetricByEggChamber('nuc',1,'Mean','raw',[6,7,8]); 

% nuclear mean intensity in DAPI channel 
% (corrected for background estimated across the egg chambers in the image)
ect.scatterPlotMetricByEggChamber('nuc',1,'Mean','wholeImgCorr','all'); 
ect.scatterPlotMetricByEggChamber('nuc',1,'Mean','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',2,'Mean','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',3,'Mean','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',4,'Mean','wholeImgCorr',[6,7,8]); 

ect.scatterPlotMetricByEggChamber('nuc',1,'Median','wholeImgCorr','all'); 
ect.scatterPlotMetricByEggChamber('nuc',1,'Median','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',2,'Median','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',3,'Median','wholeImgCorr',[6,7,8]); 
ect.scatterPlotMetricByEggChamber('nuc',4,'Median','wholeImgCorr',[6,7,8]); 

%% plot histograms

