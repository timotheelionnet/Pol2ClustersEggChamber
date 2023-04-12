
%% load data into master tables
compileSamplesIntoGlobalTables;

%% set path to include color maps used below
addpath('cbrewer');

scatter3NucleiAndClusters(nrs,cls,1,'um','all','nucVolume','clustC2_Mean_plasmCorr');


%% plot nucleus and nucloplasm intensity as a function of nucleus ID
eggChamberID = unique(nrs.eggChamberID);
ho = 0.1; % horizontal offset from figure corner for plots
vo = 0.1; % vertical offset from figure corner for plots
hs = 0.1; % horizontal spacing between plots
vs = 0.05; % vertical spacing between plots

for i=1:numel(eggChamberID)
    fNucInt = figure('Name',['egg chamber ',num2str(eggChamberID(i)),...
        ' - Nucleus and Nucleoplasm intensities']);

    [nChannels,c] = findChannelNumberFromNucTable(nrs);
    
    nRows = nChannels;
    nCols = 2;
    plotSize = [(1-2*ho-(nCols-1)*hs)/nCols,...
        (1-2*vo-(nRows-1)*vs)/nRows];
    
    for j=1:nChannels
        colIdx = 1;
        pos = [ho+(colIdx-1)*(plotSize(1)+hs),vo+(nChannels - j)*(plotSize(2)+vs),...
            plotSize(1),plotSize(2)];
        aNucInt{j} = axes(fNucInt,'Position',pos);
        if j==1
            title('Nuclear Intensity');
        end
        ylabel(['C',num2str(j),' Int']);
        hold on;

        colIdx = 2;
        pos = [ho+(colIdx-1)*(plotSize(1)+hs),vo+(nChannels - j)*(plotSize(2)+vs),...
            plotSize(1),plotSize(2)];
        aPlasmInt{j} = axes(fNucInt,'Position',pos);
        if j==1
            title('Nucleoplasm Intensity');
        end
        ylabel(['C',num2str(j),' Int']);
        hold on;
    end

    % select rows of nucleus table that correspond to current
    % egg chamber
    idxEn = nrs.eggChamberID == eggChamberID(i);
    curNrs = nrs(idxEn,:);
    
    % set color map
    nucIDs = unique(curNrs.nucID);
    nNuclei = numel(nucIDs);
    nucleiCmap = cbrewer('qual','Accent',nNuclei);

    for j=1:nChannels
        for k=1:nNuclei
            errorbar(aNucInt{j},...
                curNrs.nucID(k),...
                curNrs.(['nucC',num2str(j),'_Mean_eggChamberCorr'])(k),...
                curNrs.(['nucC',num2str(j),'_StdDev_raw'])(k),'o',...
                'Color',nucleiCmap(k,:),'MarkerFaceColor','auto',...
                'MarkerEdgeColor','auto',...
                'DisplayName',['C',num2str(j),'nuc',num2str(k),'mean/StdDev Int']);
            errorbar(aPlasmInt{j},...
                curNrs.nucID(k),...
                curNrs.(['plasmC',num2str(j),'_Mean_eggChamberCorr'])(k),...
                curNrs.(['plasmC',num2str(j),'_StdDev_raw'])(k),'o',...
                'Color',nucleiCmap(k,:),'MarkerFaceColor','auto',...
                'MarkerEdgeColor','auto',...
                'DisplayName',['C',num2str(j),'plasm',num2str(k),'mean/StdDev Int']);
        end
    end

end

%% plot cluster avg intensity and number of clusters for each nucleus
eggChamberID = unique(nrs.eggChamberID);
ho = 0.1; % horizontal offset from figure corner for plots
vo = 0.1; % vertical offset from figure corner for plots
hs = 0.1; % horizontal spacing between plots
vs = 0.05; % vertical spacing between plots

for i=1:numel(eggChamberID)
    fClustInt = figure('Name',['egg chamber ',num2str(eggChamberID(i)),...
        ' - Cluster Number and Intensities per nucleus']);

    [nChannels,c] = findChannelNumberFromNucTable(nrs);
    
    nRows = nChannels;
    nCols = 2;
    plotSize = [(1-2*ho-(nCols-1)*hs)/nCols,...
        (1-2*vo-(nRows-1)*vs)/nRows];
    
    for j=1:nChannels
        colIdx = 1;
        pos = [ho+(colIdx-1)*(plotSize(1)+hs),vo+(nChannels - j)*(plotSize(2)+vs),...
            plotSize(1),plotSize(2)];
        aClustInt{j} = axes(fClustInt,'Position',pos);
        if j==1
            title('Avg Cluster Intensity');
        end
        ylabel(['C',num2str(j),' Int']);
        hold on;

        colIdx = 2;
        if j==1
            pos = [ho+(colIdx-1)*(plotSize(1)+hs),vo+(nChannels - j)*(plotSize(2)+vs),...
                plotSize(1),plotSize(2)];
            aClustNum{j} = axes(fClustInt,'Position',pos);
            
            title('Cluster Number');
        
            ylabel('number');
            hold on;
        end
    end

    % select rows of nucleus table that correspond to current
    % egg chamber
    idxEn = nrs.eggChamberID == eggChamberID(i);
    curNrs = nrs(idxEn,:);
    
    % set color map
    nucIDs = unique(curNrs.nucID);
    nNuclei = numel(nucIDs);
    nucleiCmap = cbrewer('qual','Accent',nNuclei);

    for j=1:nChannels
        for k=1:nNuclei
            if curNrs.nClusters(k) ~=0
                errorbar(aClustInt{j},...
                    curNrs.nucID(k),...
                    curNrs.(['clustC',num2str(j),'_MeanMean_plasmCorr'])(k),...
                    curNrs.(['clustC',num2str(j),'_StdDevMean_plasmCorr'])(k),'o',...
                    'Color',nucleiCmap(k,:),'MarkerFaceColor','auto',...
                    'MarkerEdgeColor','auto',...
                    'DisplayName',['C',num2str(j),'clust',num2str(k),'mean/StdDev Int']);
            end
            if j == 1
                plot(aClustNum{j},...
                    curNrs.nucID(k),...
                    curNrs.nClusters(k),'o',...
                    'Color',nucleiCmap(k,:),'MarkerFaceColor',nucleiCmap(k,:),...
                    'MarkerEdgeColor',nucleiCmap(k,:),...
                    'DisplayName','nClusters');
                
            end
        end
    end

end



