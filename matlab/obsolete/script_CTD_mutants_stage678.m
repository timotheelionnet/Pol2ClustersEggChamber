
fijiOutFolder = '/Users/lionnt01/Dropbox/data/feiyue/20230415_CTDmutants/out';

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

% all stages
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'Volume','','nucVolume',{},[5,1,2,3,4,6],'all',0.3);

% only 7,8
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'Volume','','nucVolume',{},[5,1,2,3,4,6],[7,8],0.3);

% only 7,8, excluding 42con sample 2, EC3
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'Volume','','nucVolume',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);

%% scatter plot nuclei intensity by egg chamber: 
% raw:
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'Mean','raw','meanNucIntDAPI_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',2,'Mean','raw','meanNucIntFLAG_variantCTD_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',3,'Mean','raw','meanNucIntCycT_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',4,'Mean','raw','meanNucIntHA_refCTD_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);

% background subtracted (eggchamber):
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'Mean','eggChamberSubtracted','meanNucIntDAPI_ecCorr',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',2,'Mean','eggChamberSubtracted','meanNucIntFLAG_variantCTD_ecCorr',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',3,'Mean','eggChamberSubtracted','meanNucIntCycT_ecCorr',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',4,'Mean','eggChamberSubtracted','meanNucIntHA_refCTD_ecCorr',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);

%% scatter plot nucleoplasm intensity by egg chamber (raw): 

% nucleoplasm (raw) :
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',1,'Mean','raw','meanPlasmIntDAPI_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',2,'Mean','raw','meanPlasmIntFLAG_variantCTD_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',3,'Mean','raw','meanPlasmIntCycT_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',4,'Mean','raw','meanPlasmIntHA_refCTD_raw',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);

% nucleoplasm corrected for background estimated across the nucleolus:
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',1,'Mean','nucleoliSubtracted','meanPlasmIntDAPI_nucleoliSub',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',2,'Mean','nucleoliSubtracted','meanPlasmIntFLAG_variantCTD_nucleoliSub',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',3,'Mean','nucleoliSubtracted','meanPlasmIntCycT_nucleoliSub',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'plasm',4,'Mean','nucleoliSubtracted','meanPlasmIntHA_refCTD_nucleoliSub',{[3,2,3]},[5,1,2,3,4,6],[7,8],0.3); 

%% cluster intensity as a function of cluster volume, raw
figure('Name','Cluster Volume vs DAPI intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C1Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C1Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs FLAG (CTD variant) intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C2Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C2Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs CycT intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C3Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C3Mean plasmCorr');
set(gca,'xscale','log');

figure('Name','Cluster Volume vs HA (WT CTD) intensity');
s = scatter(ec.clustT.clust_Volume,ec.clustT.clust_C4Mean_plasmCorr,'o',...
    'MarkerFaceColor',[0.25,0.25,0.25],'MarkerEdgeColor',[0.25,0.25,0.25]);
alpha(s,0.05);
xlabel('Volume');
ylabel('cluster C4Mean plasmCorr');
set(gca,'xscale','log');

%% scatter plot by egg chamber: Number of clusters per nucleus, Cluster volume
minVolume = 1;
maxVolume = Inf;

[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveNucleusMetricByEggChamber( ...
    'nuc',1,'NumClustersMinVol','','nClusters',...
                {},[5,1,2,3,4,6],[7,8],0.3);

[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveClusterMetricByEggChamber( ...
    'clust',1,'Volume','','clusterVolume',...
                {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.3);
%% cluster intensity
minVolume =1;
maxVolume = Inf;
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveClusterMetricByEggChamber( ...
    'clust',1,'Median','nucleoliSubtracted','clustInt_DAPI_nucleoliSubtracted',...
                {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveClusterMetricByEggChamber( ...
    'clust',2,'Median','nucleoliSubtracted','clustInt_FLAG_variantCTD_nucleoliSubtracted',...
                {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveClusterMetricByEggChamber( ...
    'clust',3,'Median','nucleoliSubtracted','clustInt_CycT_nucleoliSubtracted',...
                {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.3);
[nucTable,avgNucTable,fh] = ec.scatterPlotAndSaveClusterMetricByEggChamber( ...
    'clust',4,'Median','nucleoliSubtracted','clustInt_HA_refCTD_nucleoliSubtracted',...
                {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.3);


%% same for SMALL clusters
minVolume = 0;
maxVolume = 1;
[clustTable2,avgClustTable2] = ec.scatterPlotAndSaveClusterMetricByEggChamber(...
    'clust',2,'Median','nucleoliSubtracted',...
    'clusterInt_FLAG_nucleoliSubtracted',...
    {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.1);

[clustTable4,avgClustTable4] = ec.scatterPlotAndSaveClusterMetricByEggChamber(...
    'clust',4,'Median','nucleoliSubtracted',...
    'clusterInt_HA_nucleoliSubtracted',...
    {},[5,1,2,3,4,6],[7,8],minVolume,maxVolume,0.1);

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



%% Cluster intensity, mutant vs WT allele
minVolume = 1;
idx0 = ismember(ec.clustT.eggChamber_Stage , [7,8]) ...
    & ec.clustT.clust_Volume >= minVolume;

condColors = cbrewer('qual','Set1',numel(ec.condIndices));
idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

figure('Name','FLAG (y) vs HA (x)'); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.clust_C4Median_nucleoliSubtracted(idx{i}),ec.clustT.clust_C2Median_nucleoliSubtracted(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 

alpha(0.3);
xlabel('C4Median nucleoliSubtracted');
ylabel('C2Median nucleoliSubtracted');
grid on;
legend show;

%% Cluster intensity, mutant vs WT allele, relative to nuclei levels
minVolume = 1;
idx0 = ismember(ec.clustT.eggChamber_Stage , [7,8]) ...
    & ec.clustT.clust_Volume >= minVolume;

condColors = cbrewer('qual','Set1',numel(ec.condIndices));
idx = cell(numel(ec.condIndices),1);
for i=1:numel(ec.condIndices)
    idx{i} = idx0 & ec.clustT.cond_Idx == ec.condIndices(i);
end

figure('Name','FLAG (y) vs HA (x), relative'); hold;
for i=1:numel(ec.condIndices)
    scatter(ec.clustT.clust_C4Median_nucleoliSubtracted(idx{i})./ec.clustT.plasm_C4Median_nucleoliSubtracted(idx{i}),...
        ec.clustT.clust_C2Median_nucleoliSubtracted(idx{i})./ec.clustT.plasm_C2Median_nucleoliSubtracted(idx{i}),...
        'o','MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),'DisplayName',ec.conditionNames{i});
end 
alpha(0.3);
xlabel('C4 clust Median nucleoliSubtracted, relative to plasm');
ylabel('C2 clust Median nucleoliSubtracted, relativce to plasm');
grid on;
legend show;



%% flag/ha, by sample - Normalized to nuclear levels (nucleoli Subtracted)
stages = [7,8];
condColors = cbrewer('qual','Set1',numel(ec.condIndices));
idx = ec.clustT.clust_Volume >= minVolume;
figure('Name','cluster FLAG/HA (relative to nucleoplasm levels, nucleoliSub) in clusters'); 
hold;
c = unique(ec.clustT.cond_Idx);
curX = 0;
xtl = {};
xt = [];
curCond = 0;
for i=1:numel(c)
    s = unique(ec.clustT.sample_Idx);
    for j=1:numel(s)
        e = unique(ec.clustT.eggChamber_Idx);
        for k=1:numel(e)
            curIdx = idx & ec.clustT.cond_Idx == c(i) & ec.clustT.sample_Idx == s(j) ...
                & ec.clustT.eggChamber_Idx == e(k) ...
                & ismember(ec.clustT.eggChamber_Stage,stages);
            if sum(curIdx)>0
                 flag = ec.clustT.clust_C2Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C2Median_nucleoliSubtracted(curIdx);        
                 ha = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdx);  
                 scatter(curX*ones(size(flag))-0.15+0.3*rand(numel(flag),1),flag./ha,'o',...
                        'MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),...
                        'DisplayName',[ec.conditionNames{i},' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 xt(end+1) = curX;
                 xtl{end+1} = [ec.conditionNames{i},' smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 curX = curX+1;
            end
        end
    end
end
alpha(0.3);
ylabel('flag/ha');
xticks(xt);
xticklabels(xtl);
xtickangle(45);


%% flag/ha, by sample - SMALL CLUSTERS!!! Normalized to nuclear levels (nucleoli Subtracted)
stages = [7,8];
condColors = cbrewer('qual','Set1',numel(ec.condIndices));
idx = ec.clustT.clust_Volume < minVolume;
figure('Name','cluster FLAG/HA (relative to nucleoplasm levels, nucleoliSub) in SMALL clusters'); 
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
                 flag = ec.clustT.clust_C2Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C2Median_nucleoliSubtracted(curIdx);        
                 ha = ec.clustT.clust_C4Median_nucleoliSubtracted(curIdx)./ec.clustT.plasm_C4Median_nucleoliSubtracted(curIdx);  
                 scatter(curX*ones(size(flag))-0.15+0.3*rand(numel(flag),1),flag./ha,'o',...
                        'MarkerFaceColor',condColors(i,:),'MarkerEdgeColor',condColors(i,:),...
                        'DisplayName',[ec.conditionNames{i},' Sample ',num2str(s(j)),' egg chamber ',num2str(e(k))]);
                 xt(end+1) = curX;
                 xtl{end+1} = [ec.conditionNames{i},' smpl ',num2str(s(j)),' eggChbr ',num2str(e(k))];
                 curX = curX+1;
            end
        end
    end
end
alpha(0.3);
ylabel('flag/ha');
xticks(xt);
xticklabels(xtl);
xtickangle(45);


%% PolII, by sample, in !!!!!!SMALL!!!!!  MPM2+/- clusters - Normalized to nuclear levels (nucleoli Subtracted)
minMPM2 = 2500;
stages = 10;

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
stages = 10;

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
stages = 10;

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
stages = 10;

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
stages = 10;

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
stages = 10;

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
stages = 10;

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
stages = 10;

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
