fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
fijiOutFolder = '/Users/lionnt01/Documents/data/feiyue/clusterImgAnalysis/out';
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/out_v2';

%% set path for subfunctions

addpath('subfunctions/');
addpath('cbrewer/');
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

%%
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Mean','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Mean','raw',[6,7,8]);

ec.scatterPlotNucleiMetricByEggChamber('nucleoli',1,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',2,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',3,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nucleoli',4,'Median','raw',[6,7,8]);

%% scatter plot by egg chamber: Nucleus Intensity
% nuclear mean intensity in all channels

% raw:
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','raw','all'); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Median','raw',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Median','raw',[6,7,8]);

% background subtracted (whole image)
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','wholeImgSubtracted','all'); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'Median','wholeImgSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',2,'Median','wholeImgSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',3,'Median','wholeImgSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('nuc',4,'Median','wholeImgSubtracted',[6,7,8]); 

% background subtracted (sampleROI)
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Median','sampleROISubtracted','all'); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Median','sampleROISubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Median','sampleROISubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Median','sampleROISubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Median','sampleROISubtracted',[6,7,8]); 

% nucleoplasm corrected for background estimated across the nucleolus:
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Median','nucleoliSubtracted','all'); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',1,'Median','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Median','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Median','nucleoliSubtracted',[6,7,8]); 
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Median','nucleoliSubtracted',[6,7,8]); 

%% 
ec.scatterPlotNucleiMetricByEggChamber('nucAvgClustMinVol',2,'Mean','raw','all'); 

%% 
yMetric = 'nucAvgClustMinVol_C3Mean_nucleoliSubtracted_cT1';
fh = ec.scatterPlotArbitraryMetricByEggChamber(0,'nuc',yMetric,idx,eggChamberStagesToInclude)

%%
minVolume =1;
ec.scatterPlotNucleiMetricByEggChamber('plasm',2,'Median','raw',[6,7,8]);
ec.scatterPlotNucleiMetricByEggChamber('plasm',3,'Median','raw',[6,7,8]);
ec.scatterPlotNucleiMetricByEggChamber('plasm',4,'Median','raw',[6,7,8]);

ec.scatterPlotClusterMetricByEggChamber('clust',2,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',3,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',4,'Median','raw',[6,7,8],minVolume);

%% scatter plot by egg chamber: Number of clusters per nucleus

ec.scatterPlotNucleiMetricByEggChamber('nuc',1,'NumClustersMinVol','',[6,7,8]); 

ec.scatterPlotClusterMetricByEggChamber('clust',1,'Volume','',[6,7,8],minVolume);
%%
figure('Name','Volume vs DAPI intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C1Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('C1Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Volume vs MPM2 intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C2Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('C2Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Volume vs Ser5ph intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C3Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('C3Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Volume vs Pol II intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C4Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('C4Mean plasmCorr');
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
ec.scatterPlotClusterMetricByEggChamber('clust',1,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',2,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',3,'Median','raw',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',4,'Median','raw',[6,7,8],minVolume);

ec.scatterPlotClusterMetricByEggChamber('clust',1,'Median','plasmCorr',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',2,'Median','plasmCorr',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',3,'Median','plasmCorr',[6,7,8],minVolume); 
ec.scatterPlotClusterMetricByEggChamber('clust',4,'Median','plasmCorr',[6,7,8],minVolume); 

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

%% MPM2/PolII, by sample, in MPM2+ clusters / Normalized to nuclear levels (wholeImgSubtracted)
minMPM2 = 5000;
stages = [6,7,8];

idx = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','S5ph/Pol2 (relative to nuclear levels wholeImgSub) in MPM2+ clusters'); 
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
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_wholeImgSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_wholeImgSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
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
figure('Name','S5ph/Pol2 (relative to nuclear levels wholeImgSub) in MPM2- clusters'); 
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
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_wholeImgSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_wholeImgSubtracted(curIdx);  
                 if i==1
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0.5,0,0],'MarkerEdgeColor',[0.5,0,0],...
                        'DisplayName',['Cond ',num2str(c(i)),' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 else
                    scatter(curX*ones(size(s5))-0.15+0.3*rand(numel(s5),1),s5./p2,'o','MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5],...
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

%% MPM2/PolII, by sample, in MPM2+ clusters / Normalized to nuclear levels (sampleROISubtracted)
minMPM2 = 2500; 
stages = [6,7,8];

idx = ec.clustT.clust_C2Median_plasmCorr > minMPM2...
    & ec.clustT.clust_Volume >= minVolume;
figure('Name','S5ph/Pol2 (relative to nuclear levels sampleROI) in MPM2+ clusters'); 
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
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_sampleROISubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_sampleROISubtracted(curIdx);  
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
figure('Name','S5ph/Pol2 (relative to nuclear levels sampleROI) in MPM2- clusters'); 
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
                 s5 = ec.clustT.clust_C3Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C3Median_wholeImgSubtracted(curIdx);        
                 p2 = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_wholeImgSubtracted(curIdx);  
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


%% MPM2/PolII, by sample, in MPM2+ clusters / NOT Normalized to nuclear levels
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
