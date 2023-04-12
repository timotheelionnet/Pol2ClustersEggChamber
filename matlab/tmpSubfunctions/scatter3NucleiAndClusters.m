function scatter3NucleiAndClusters(nrs,cls,eggChamberID,DistUnits,nucIDs,nucColorVar,clustColorVar)
% plots in 3D the position of the nuclei and associated clusters for a single egg chamber.
% Color code can be whatever variable.

% input:
% nrs: simplified nucleus metrics table
% cls: simplified cluster metrics table
% eggChamberID: index of the desired eggChamber 
    % (needs to match a value of the eggChamberID in the metrics tables)
% distUnits: char used to write the distance units on the axes (e.g. 'um')
% nucIDs: nucID values of the nuclei to be plotted. 
    % enter either 'all' (all nuclei in the egg chamber), or a single value
    % (e.g. 1), or an array of several values ([1,7,8])
% nucColorVars: variable to be used to color-code  
    % (needs to be a variable in the nrs table)
% clustColorVars: variable to be used to color-code clusters 
    % (needs to be a variable in the cls table)

% TO DO
% add choice color map assigned to nuclei or cluster cData
% get the colormap to plot correctly and assign ticklabels
% find a way to display one or the others in a neutral color.
% make sure the case where the colors are single-valued doesnt give an
% error

% build figure
fh = figure('Name',['egg chamber ',num2str(eggChamberID),...
        ' - Positions of nuclei and clusters'],'Position',[300,500,1000,500],'Units','points');
aC = axes(fh,'Position',[0.17,0.1,0.65,0.5]);
aC.Visible = 'off';
aN = axes(fh,'Position',[0.17,0.1,0.65,0.85]);
hold;

% select rows of nucleus and cluster tables that correspond to selected
% egg chamber
idxE = nrs.eggChamberID == eggChamberID;
curNrs = nrs(idxE,:);
idxE = cls.eggChamberID == eggChamberID;
curCls = cls(idxE,:);
maxNuclei = max(curNrs.nucID);

% select rows of nucleus and cluster tables that correspond to selected 
% nucleus
%construct list of nucIDs if 'all' option selected
if ischar(nucIDs)
    if strcmp(nucIDs,'all')
        nucIDs = unique(curNrs.nucID);
    else
        disp('could not parse nucIDs entry');
    end
end

% keep rows of the table that have the eggchamber and nuclei of interest
idxN = ismember(curNrs.nucID,nucIDs);
curNrs = curNrs(idxN,:);
idxN = ismember(curCls.nucID,nucIDs);
curCls = curCls(idxN,:);

%% plot nuclei

% build nucleiCmap, a nNuclei x 3 array that lists the RGB color code triplet for each
% nucleus in the list
nNuclei = numel(nucIDs);
if strcmp(nucColorVar,'nucID')
    nucColors = cbrewer('qual','Accent',maxNuclei);
    nucCmap = [];
    for i=1:size(curNrs,1)
        nucCmap = [nucCmap;nucColors(curNrs.nucID(i),: )];
    end
    colormap(aN,nucColors);
    % set up nuclei colorbar
    nucTickLabels = cellfun(@num2str,num2cell(1:maxNuclei),'UniformOutput',0);
    cN = colorbar(aN,'Ticks',((1:maxNuclei)-0.5)/maxNuclei,...
        'Ticklabels',nucTickLabels,...
        'Position',[.07 .1 .03 .8]);
    cN.Label.String = ['nucleus color: ',strrep(nucColorVar,'_','\_')];

else
    nucColors = cbrewer('div','RdYlBu',255);
    if max(curNrs.(nucColorVar)) == min(curNrs.(nucColorVar))
        colorVals = 123*ones(size(curNrs.(nucColorVar)));
    else
        colorVals = 1+ceil(254* ...
            (curNrs.(nucColorVar)-min(curNrs.(nucColorVar))) ...
            /(max(curNrs.(nucColorVar)) - min(curNrs.(nucColorVar))));
    end
    nucCmap = nucColors(colorVals,:);

    % set up nuclei colorbar
    colormap(aN,nucColors);
    cN = colorbar(aN,'Position',[.07 .1 .03 .8]);
    if max(curNrs.(nucColorVar)) ~= min(curNrs.(nucColorVar))
        clim([min(curNrs.(nucColorVar)),max(curNrs.(nucColorVar))]);
    elseif max(curNrs.(nucColorVar)) > 0
        clim([0,2*max(curNrs.(nucColorVar))]);
    elseif max(curNrs.(nucColorVar)) == 0
        clim([-1,1]);
    else
        max(curNrs.(nucColorVar))
        clim([2*max(curNrs.(nucColorVar)),0]);
    end
    cN.Label.String = ['nucleus color: ',strrep(nucColorVar,'_','\_')];
end

% plot nuclei
scatter3(aN,curNrs.nucCentroid_X,curNrs.nucCentroid_Y,curNrs.nucCentroid_Z,...
    5000,nucCmap,'filled','o');

xlabel(['X in ',DistUnits]);
ylabel(['Y in ',DistUnits]);
zlabel(['Z in ',DistUnits]);
grid on

%% plot clusters

% build clustCmap, a Nclusters x 3 arrays that lists the RGB color code triplet for
% each cluster
if strcmp(clustColorVar,'nucID')
    nucColors = cbrewer('qual','Accent',maxNuclei);
    clustCmap = [];
    for i=1:size(curCls,1)
        clustCmap = [clustCmap;nucColors(curCls.nucID(i),: )];
    end

    % set up clusters colorbar
    clustTickLabels = cellfun(@num2str,num2cell(1:maxNuclei),'UniformOutput',0);
    colormap(aC,nucColors);
    cC = colorbar(aC,'Ticks',((1:maxNuclei)-0.5)/maxNuclei,...
        'Ticklabels',clustTickLabels,...
        'Position',[.9 .1 .03 .8]);
    cC.Label.String = ['clusters color: ',strrep(clustColorVar,'_','\_')];
else
    clustColors = cbrewer('div','RdYlBu',255);
    if max(curCls.(clustColorVar)) == min(curCls.(clustColorVar))
        colorVals = 123*ones(size(curCls.(clustColorVar)));
    else
        colorVals = 1+ceil(254*...
            (curCls.(clustColorVar) - min(curCls.(clustColorVar)))...
            /(max(curCls.(clustColorVar)) - min(curCls.(clustColorVar))));
    end
    clustCmap = clustColors(colorVals,:);

    % set up clusters colorbar
    colormap(aC,clustColors);
    cC = colorbar(aC,'Position',[.9 .1 .03 .8]);
    axes(aC);
    if max(curCls.(clustColorVar)) ~= min(curCls.(clustColorVar))
        clim([min(curCls.(clustColorVar)),max(curCls.(clustColorVar))]);
    elseif max(curCls.(clustColorVar)) > 0
        clim([0,2*max(curCls.(clustColorVar))]);
    elseif max(curCls.(clustColorVar)) == 0
        clim([-1,1]);
    else
        clim([2*max(curCls.(clustColorVar)),0]);
    end
    cC.Label.String = ['clusters color: ',strrep(clustColorVar,'_','\_')];
end
axes(aN);
scatter3(aN,curCls.clustCentroid_X,curCls.clustCentroid_Y,curCls.clustCentroid_Z,...
        100,clustCmap,'filled','o');
alpha 0.5
fontsize(fh,14,'points');


