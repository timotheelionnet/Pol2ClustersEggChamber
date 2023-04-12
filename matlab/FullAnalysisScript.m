fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/testSer5ph/nucSeg20_3img';

% intitialize object 
ec = eggChamberDataFolder(fijiOutFolder);

% load raw data from nuclei segmentation
t = ec.loadAllEggChamberData();

%% compute background-corrected nuclei intensity values (eggChamber and whole image)
t = backgroundCorrectNucIntensity(t);

%% assign nuclei to groups
t = assignNucleiToGroupsWholeTable(t);

%% plot nuclei stats by group
plotNucleiStatsByGroup(t,'nucVolume');

t = addvars(t,t.nucC1_Mean_wholeImgCorr.*t.nucVolume,'NewVariableNames','integratedDAPI');
plotNucleiStatsByGroup(t,'integratedDAPI');
%%