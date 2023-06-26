fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/out_v4';

%% set path for subfunctions

addpath('subfunctions/');
addpath('cbrewer/');
%%

% intitialize eggChamberDataFolder objec - this collecs the locations of
% the files associated with the experiment and stores them by condition and sample
ec = eggChamberDataSet(fijiOutFolder);

%% load nucleus-wide raw data
% load raw data from nuclei segmentation - builds a table e
% c.nucFullT that combines
% nuclei from all samples (one row per nucleus) and compiles dozens of
% metrics of size and intensity in every channel. (creates nucT as a backup copy)
ec.loadAllEggChamberNucleiData;

%% load cluster data
% populates the table ec.clustT
% also adds nucleoli and nucleoplasm data to ec.nucT and ec.nucFullT tables
% does NOT add summary cluster metrics to nucT and nucFullT
ec.loadAllClusterData();

%% background correct nuclei intensity
% subtracts image-wide background value (wholeImgSubtracted) and typical background
% in an eggchamber away from nuclei (sampleROISubtracted)
ec.backgroundCorrectNucIntensity;

%% add nuc stats to cluster table
% takes metrics from nucFullT and copies them as extra variables of the
% cluster table.
%ec.clustT 
ec.clustT = ec.addNucStatsToClustTable();

%% background correct clusters intensity
% subtracts nucleoplasm intensity (plasmSubtracted) and nucleoli intensity
% (nucleoliSubtracted)
ec.backgroundCorrectClustIntensity();

%% add cluster data to nuc table
% takes metrics from clustT, averages them over each nucleus
% and copies them as extra variables of the nuc table.
ec.addAverageClusterStatsToNucTable();

%% scatter plot by egg chamber: Nucleus Volume
% plot where all nuclei from an egg chamber are grouped separately.
ec.scatterPlotNucleiMetricByEggChamber('eggChamber',1,'Stage','','all'); % egg chmaber stage
ec.scatterPlotNucleiMetricByEggChamber('eggChamber',1,'Stage','',[6,7,8]); % egg chmaber stage
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Volume','','all'); % nuclear Volume
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Volume','',[6,7,8]); % nuclear Volume
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Volume','',[6,7]); % nuclear Volume


%% scatter plot nuclei intensity by egg chamber: 
% raw:
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Mean','raw',[6,7,8]);

% displaying table in command window with the raw values estimated for four
% nuclei in ctrl_1002 image to match w visual estimates.
% variables in order are nucleoli (raw),nucleoplasm (raw) for all four
% channels.
cond = 1;
sample = 2;
nucIdx = [20,3,9,7];
b = zeros(numel(nucIdx),4);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    
    b(i,:) = a;
end
disp('raw nuclei intensities for four nuclei in Ctrl 1002 FOV:');
b

% background subtracted (eggchamber)
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Mean','eggChamberSubtracted',[6,7,8]);

% displaying table in command window with the raw values estimated for four
% nuclei in ctrl_1002 image to match w visual estimates.
% variables in order are nucleoli (raw),nucleoplasm (raw) for all four
% channels.
cond = 1;
sample = 2;
nucIdx = [20,3,9,7];
b = zeros(numel(nucIdx),4);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['nuc_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    
    b(i,:) = a;
end
disp('eggChamber-subtracted nuclei intensities for four nuclei in Ctrl 1002 FOV:');
b

%% scatter plot nucleoplasm intensity by egg chamber (raw): 

% nucleoplasm (raw) :
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Mean','raw',[6,7,8]); 

% nucleoli (raw) :
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',1,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',2,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',3,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',4,'Mean','raw',[6,7,8]); 

% displaying table in command window with the raw values estimated for four
% nuclei in ctrl_1002 image to match w visual estimates.
% variables in order are nucleoli (raw),nucleoplasm (raw) for all four
% channels.
cond = 1;
sample = 2;
nucIdx = [20,3,9,7];
b = zeros(numel(nucIdx),8);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    
    b(i,:) = a;

end
b

%%
% nucleoplasm corrected for background estimated across the nucleolus:
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Mean','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Mean','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Mean','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Mean','nucleoliSubtracted',[6,7,8]); 

% displaying table in command window with the values estimated for four
% nuclei in ctrl_1002 image to match w visual estimates.
% variables in order are nucleoplasm (nucleoli corrected) for all four
% channels.
cond = 1;
sample = 2;
nucIdx = [20,3,9,7];
b = zeros(numel(nucIdx),4);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_nucleoliSubtracted'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_nucleoliSubtracted'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_nucleoliSubtracted'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_nucleoliSubtracted'])(imgIdx)];
    
    b(i,:) = a;
end
disp('nucleoli-subtracted nucleoplasm intensities for four nuclei in Ctrl 1002 FOV:');
b

% nucleoplasm corrected for background estimated across each egg chamber
% outside of nuclei
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Mean','eggChamberSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Mean','eggChamberSubtracted',[6,7,8]); 

% displaying table in command window with the raw values estimated for four
% nuclei in ctrl_1002 image to match w visual estimates.
% variables in order are nucleoplasm (eggChmaber-corrected) for all four
% channels.
cond = 1;
sample = 2;
nucIdx = [20,3,9,7];
b = zeros(numel(nucIdx),4);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_eggChamberSubtracted'])(imgIdx)];
    
    b(i,:) = a;
end
disp('eggChamber-subtracted nucleoplasm intensities for four nuclei in Ctrl 1002 FOV:');
b


%%
minVolume =1;
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Median','raw',[6,7,8]);
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Median','raw',[6,7,8]);
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Median','raw',[6,7,8]);

ec.scatterPlotClusterMetricByEggChamber('clust',2,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',3,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',4,'Median','raw',[6,7,8],minVolume);

%% scatter plot by egg chamber: Number of clusters per nucleus
minVolume =1;
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'NumClustersMinVol','',[6,7,8]); 

ec.scatterPlotClusterMetricByEggChamber('clust',1,'Volume','',[6,7,8],minVolume);

%% 
ec.scatterPlotNucleiMetricByEggChamber('nucAvgClustMinVol',2,'Mean','raw','all'); 

%% 
yMetric = 'nucAvgClustMinVol_C3Mean_nucleoliSubtracted_cT1';
fh = ec.scatterPlotArbitraryMetricByEggChamber(0,'nuc',yMetric,idx,eggChamberStagesToInclude)

%% cluster intensity as a function of cluster volume, raw
figure('Name','Cluster Volume vs DAPI intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C1Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C1Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs MPM2 intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C2Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C2Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs Ser5ph intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C3Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C3Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs Pol II intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C4Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C4Mean plasmCorr');
set(gca,'xscale','log');

%% cluster intensity as a function of cluster volume, normalized by nucleoplasm intensity
figure('Name','Cluster Volume vs log2 DAPI intensity');
s = scatter(ec.clustT.clust_Volume,log(ec.clustT.clust_C1Mean_plasmCorr./ec.clustT.plasm_C1Mean_eggChamberSubtracted)/log(2),'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('log2 (cluster C1Mean plasmCorr / plasm C1Mean eggChamberSubtracted)');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs log2 MPM2 intensity');
s = scatter(ec.clustT.clust_Volume,log(ec.clustT.clust_C2Mean_plasmCorr./ec.clustT.plasm_C2Mean_eggChamberSubtracted)/log(2),'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('log2 (cluster C2Mean plasmCorr / plasm C2Mean  eggChamberSubtracted)');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs log2 Ser5ph intensity');
s = scatter(ec.clustT.clust_Volume,log(ec.clustT.clust_C3Mean_plasmCorr./ec.clustT.plasm_C3Mean_eggChamberSubtracted)/log(2),'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('log2 (cluster C3Mean plasmCorr / plasm C3Mean  eggChamberSubtracted)');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs log2 Pol II intensity');
s = scatter(ec.clustT.clust_Volume,log(ec.clustT.clust_C4Mean_plasmCorr./ec.clustT.plasm_C4Mean_eggChamberSubtracted)/log(2),'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('log2 (cluster C4Mean plasmCorr / / plasm C4Mean  eggChamberSubtracted)');
set(gca,'xscale','log');
%% scatter plot clusters Volume
minVolume = 1;

% all clusters
ec.scatterPlotClusterMetricByEggChamber('clust',1,'Volume','',[6,7,8],0);

% only bonafide HLBs (> 1um^3)
ec.scatterPlotClusterMetricByEggChamber('clust',1,'Volume','',[6,7,8],minVolume);

% all clusters
ec.scatterPlotClusterMetricByEggChamber('clust',1,'Mean','plasmCorr',[6,7,8],0);

% only bonafide HLBs (> 1um^3)
ec.scatterPlotClusterMetricByEggChamber('clust',1,'Mean','plasmCorr',[6,7,8],minVolume);

%% scatter plot clusters intensity

ec.scatterPlotClusterMetricByEggChamber('clust',1,'Median','nucleoliSubtracted',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',2,'Median','nucleoliSubtracted',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',3,'Median','nucleoliSubtracted',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',4,'Median','nucleoliSubtracted',[6,7,8],minVolume); 

%% plasmCorr vs plamsMean levels
minVolume = 1;
idx = ismember(ec.clustT.eggChamber_Stage , [6,7,8]) ...
    & ec.clustT.clust_Volume >= minVolume;

idx1 = idx & ec.clustT.cond_Idx == 1;
idx2 = idx & ec.clustT.cond_Idx == 2;

figure('Name','plasmCorr vs plasm_Mean, C4'); hold;
scatter(ec.clustT.plasm_C4Mean_raw(idx1),ec.clustT.clust_C4Mean_raw(idx1) - ec.clustT.clust_C4Mean_plasmCorr(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.plasm_C4Mean_raw(idx2),ec.clustT.clust_C4Mean_raw(idx2) - ec.clustT.clust_C4Mean_plasmCorr(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('plasmCorr');
ylabel('plasm Mean');
grid on;
legend show;

figure('Name','plasmCorr vs plasm_Mean, C3'); hold;
scatter(ec.clustT.plasm_C3Mean_raw(idx1),ec.clustT.clust_C3Mean_raw(idx1) - ec.clustT.clust_C3Mean_plasmCorr(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.plasm_C3Mean_raw(idx2),ec.clustT.clust_C3Mean_raw(idx2) - ec.clustT.clust_C3Mean_plasmCorr(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('plasmCorr');
ylabel('plasm Mean');
grid on;
legend show;

%% MPM2 vs Ser5ph
minVolume = 1;
idx = ismember(ec.clustT.eggChamber_Stage , [6,7,8]) ...
    & ec.clustT.clust_Volume >= minVolume;

idx1 = idx & ec.clustT.cond_Idx == 1;
idx2 = idx & ec.clustT.cond_Idx == 2;

figure('Name','MPM2 vs Ser5ph'); hold;
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx1),ec.clustT.clust_C3Median_nucleoliSubtracted(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx2),ec.clustT.clust_C3Median_nucleoliSubtracted(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('C2Median nucleoliSubtracted');
ylabel('C3Median nucleoliSubtracted');
grid on;
legend show;

figure('Name','MPM2 vs PolII'); hold;
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx1),ec.clustT.clust_C4Median_nucleoliSubtracted(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx2),ec.clustT.clust_C4Median_nucleoliSubtracted(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('C2Median nucleoliSubtracted');
ylabel('C4Median nucleoliSubtracted');
grid on;
legend show;

%% MPM2 vs Ser5ph, relative to nuclei levels
minVolume = 1;
idx = ismember(ec.clustT.eggChamber_Stage , [6,7,8]) ...
    & ec.clustT.clust_Volume >= minVolume;

idx1 = idx & ec.clustT.cond_Idx == 1;
idx2 = idx & ec.clustT.cond_Idx == 2;

figure('Name','MPM2 vs Ser5ph, relative'); hold;
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx1),ec.clustT.clust_C3Median_nucleoliSubtracted(idx1)./ec.clustT.plasm_C3Median_nucleoliSubtracted(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx2),ec.clustT.clust_C3Median_nucleoliSubtracted(idx2)./ec.clustT.plasm_C3Median_nucleoliSubtracted(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('C2Median nucleoliSubtracted');
ylabel('C3Median nucleoliSubtracted relative to nucleus');
grid on;
legend show;

figure('Name','MPM2 vs PolII, relative'); hold;
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx1),ec.clustT.clust_C4Median_nucleoliSubtracted(idx1)./ec.clustT.plasm_C4Median_nucleoliSubtracted(idx1),...
    'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],'DisplayName','Ctrl');
scatter(ec.clustT.clust_C2Median_nucleoliSubtracted(idx2),ec.clustT.clust_C4Median_nucleoliSubtracted(idx2)./ec.clustT.plasm_C4Median_nucleoliSubtracted(idx2),...
    'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],'DisplayName','TRI');
alpha(0.1);
xlabel('C2Median nucleoliSubtracted');
ylabel('C4Median nucleoliSubtracted relative to nucleus');
grid on;
legend show;

%% MPM2/PolII, by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idx = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster S5ph/Pol2 (relative to nucleoplasm levels nucleoliSub) in MPM2+ clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdx = idx & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdx)>0
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 end
                 
                 xt(end+1) = curX;
                 xtl{end+1} = ['cond ',num2str(c(i)),' smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 curX = curX+1;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph/pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

idx = ec.clustT.clust_C2Median_plasmCorr <= minMPM2 ...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster S5ph/Pol2 (relative to nucleoplasm levels nucleoliSub) in MPM2- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdx = idx & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdx)>0
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 end

                 xt(end+1) = curX;
                 xtl{end+1} = ['cond ',num2str(c(i)),' smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 curX = curX+1;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph/pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

%% PolII, by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster Pol2 (relative to nucleoplasm levels nucleoliSub) in MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 p2Pos = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxPos);  
                 if i==1
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 p2Neg = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxNeg);  
                 if i==1
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

% Ser5ph  by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster S5ph (relative to nucleoplasm levels nucleoliSub) in MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 s5Pos = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxPos);        
                 if i==1
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 s5Neg = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxNeg);        
                 if i==1
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

%% PolII, by sample, in !!!!!!SMALL!!!!!  MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume < minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume < minVolume;
figure('Name','SMALL cluster Pol2 (relative to nucleoplasm levels nucleoliSub) in MPM2+/- small clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 p2Pos = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxPos);  
                 if i==1
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 p2Neg = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxNeg);  
                 if i==1
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

% Ser5ph  by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume < minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume < minVolume;
figure('Name','SMALL cluster S5ph (relative to nucleoplasm levels nucleoliSub) in SMALL MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 s5Pos = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxPos);        
                 if i==1
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 s5Neg = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxNeg);        
                 if i==1
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

%% PolII, by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','plasm Pol2 in MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 p2Pos = ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxPos);  
                 if i==1
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 p2Neg = ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxNeg);  
                 if i==1
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

% Ser5ph  by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','plasm S5ph (nucleoliSub) in MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 s5Pos = ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxPos);        
                 if i==1
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 s5Neg = ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxNeg);        
                 if i==1
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

%% PolII, by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
% SMALL CLUSTERS!!!!
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume < minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume < minVolume;
figure('Name','plasm Pol2 in SMALL MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 p2Pos = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxPos);  
                 if i==1
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(p2Pos))-0.15+0.3*rand(numel(p2Pos),1),p2Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 p2Neg = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxNeg);  
                 if i==1
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(p2Neg))-0.15+0.3*rand(numel(p2Neg),1),p2Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);

% Ser5ph  by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','SMALL cluster S5ph (relative to nucleoplasm levels nucleoliSub) in MPM2+/- SMALL clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 s5Pos = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxPos);        
                 if i==1
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 s5Neg = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxNeg);        
                 if i==1
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph');
xticks(xt);
xticklabels(xtl);
xtickangle(45);


%% combined MPM2/PolII, by sample, in MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = [6,7,8];

idxPos = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
idxNeg = ec.clustT.clust_C2Median_plasmCorr <= minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster S5ph/Pol2 (relative to nucleoplasm levels nucleoliSub) in MPM2+/- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdxPos = idxPos & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            curIdxNeg = idxNeg & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdxPos | curIdxNeg)>0
                 s5Pos = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxPos);        
                 p2Pos = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxPos)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxPos);  
                 if i==1
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos./p2Pos,'o','MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 else
                    scatter(curX*ones(size(s5Pos))-0.15+0.3*rand(numel(s5Pos),1),s5Pos./p2Pos,'o','MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8],...
                        'DisplayName',['Ctrl Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2+']);
                 end

                 s5Neg = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C3Median_nucleoliSubtracted(curIdxNeg);        
                 p2Neg = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdxNeg)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdxNeg);  
                 if i==1
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg./p2Neg,'o','MarkerFaceColor',[0.2,0,0],'MarkerEdgeColor',[0.2,0,0],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 else
                    scatter((curX+2)*ones(size(s5Neg))-0.15+0.3*rand(numel(s5Neg),1),s5Neg./p2Neg,'o','MarkerFaceColor',[0,0,0.2],'MarkerEdgeColor',[0,0,0.2],...
                        'DisplayName',['TRI Sample ',num2str(s(j)),' egg chamber ',num2str(e(k)),' MPM2-']);
                 end
                 
                 xt(end+1) = curX+1;
                 if i==1
                    xtl{end+1} = ['Ctrl smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 else
                     xtl{end+1} = ['TRI smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 end
                 curX = curX+6;
            end
        end
    end
end
alpha(0.1);
ylabel('S5ph/pol2');
xticks(xt);
xticklabels(xtl);
xtickangle(45);



%% #OLD MPM2/PolII, by sample, in MPM2+ clusters / NOT Normalized to nuclear levels
minMPM2 = 2500; 
stages = [6,7,8];

idx = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','S5ph/Pol2 (not norm rel to nucl) in MPM2+ clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdx = idx & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdx)>0
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))+0.2*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))+0.2*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 end
                 curX = curX+1;
            end
        end
    end
end
alpha(0.1);
xlabel('Egg Chamber');
ylabel('S5ph/pol2');

idx = ec.clustT.clust_C2Median_plasmCorr <= minMPM2 ...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','S5ph/Pol2 (not norm rel to nucl) in MPM2- clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdx = idx & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdx)>0
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))+0.2*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))+0.2*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 end
                 curX = curX+1;
            end
        end
    end
end
alpha(0.1);
xlabel('Egg Chamber');
ylabel('S5ph/pol2');


%% LEFTOVERS
figure;


%% scatter plot by sample
% plot where all nuclei from each sample are grouped separately.
% not terribly useful unless the egg chambers havent been segmented yet.
% ec.scatterPlotNucleiMetricBySample('nuc',1,'Volume','');
% ec.scatterPlotNucleiMetricBySample('nuc',1,'Mean','wholeImgSubtracted');
