fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img';

% intitialize eggChamberDataFolder object - this collects the locations of
% the files associated with the experiment and stores them by condition and sample
ec = eggChamberDataFolder(fijiOutFolder);

%% load nucleus-wide raw data
% load raw data from nuclei segmentation - builds a table t that combines
% nuclei from all samples (one row per nucleus) and compiles dozens of
% metrics of size and intensity in every channel
t = ec.loadAllEggChamberData();

%% compute background-corrected nuclei intensity values (eggChamber and whole image)
t = backgroundCorrectNucIntensity(t);

%% plot raw nucleus-wide intensity stats for the different channels
% TO DO
% plotNucleiIntensityStats
% for each of raw intensity, eggchamberCorr, wholeImgCorr,
% loop through channel
% also plot size
% plot metric of cell cycle for the cycE channel

%% assign nuclei to cell cycle groups
%t = assignNucleiToGroupsWholeTable(t);

%% plot nuclei stats by group
plotNucleiStatsByGroup(t,'nucVolume');

t = addvars(t,t.nucC1_Mean_wholeImgCorr.*t.nucVolume,'NewVariableNames','integratedDAPI');
plotNucleiStatsByGroup(t,'integratedDAPI');
%%