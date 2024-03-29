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

%% compute background-corrected nuclei intensity values (eggChamber and whole image)
t2 = backgroundCorrectNucIntensity(t);

%% plot raw nucleus-wide intensity stats for the different channels
% TO DO
% plotNucleiIntensityStats
% for each of raw intensity, eggchamberCorr, wholeImgCorr,
% loop through channel
% also plot size
% plot metric of cell cycle for the cycE channel

prefix = 'nuc';
geomMetricList = {'Volume','Centroid_Z'};
channelMetricList = {'Mean'};
processingList = {'raw','wholeImgCorr','eggChamberCorr'};
spacingUnit = 1; 
freeSpaceBetweenSamples = 0.2;

for i=1:numel(geomMetricList)
    % build figure
    figure('Name',['Nuclei ',geomMetricList{i}]);
    hold;

    % generate column name
    curVarName = [prefix,geomMetricList{i}];

    % collect the x values to plot each condition/sample at
    [xSampleVals,xSampleIDs] = getSampleXValues(...
        ec.conditionNames,ec.nSamples,spacingUnit);
    ySampleVals = zeros(numel(xSampleVals),1);

    % collect the values of the metric for each condition/sample
    xPlot = [];
    yPlot = [];
    xSampleValsVec = [];
    xSampleIDsVec = {};

    ctr = 0;
    for j=1:ec.nConditions
        for k=1:ec.nSamples(j)
            ctr = ctr+1;
            % collect values for the desired metric from all nuclei for the
            % current sample/condition
            x = t.(curVarName)(t.condIdx ==j & t.sampleIdx == k);

            % generate slightly offset x coordinates for each nucleus,
            % centered around the sample X
            nNuclei = size(x,1);
            if nNuclei >1
                % spacing between nuclei
                nucSpacing = ...
                    spacingUnit*(1-2*freeSpaceBetweenSamples)/(nNuclei-1);

                % x coordinate for each nucleus of current condition/sample
                curXPlot = xSampleVals(j,k) - floor(nNuclei/2)*nucSpacing ...
                    + (0:(nNuclei-1))*nucSpacing;

                % y coordinate for each nucleus of current condition/sample
                curYPlot = t.(curVarName)(t.condIdx ==j & t.sampleIdx == k)';

            elseif nNuclei == 1
                curXPlot = xSampleVals(i,j);
                curYPlot = t.(curVarName)(t.condIdx ==j & t.sampleIdx == k)';

            elseif nNuclei == 0
                curXPlot = [];
                curYPlot = [];
            end

            % append coordinates of current condition/sample to global list
            xPlot = [xPlot,curXPlot];
            yPlot = [yPlot,curYPlot];

            xSampleValsVec = [xSampleValsVec,xSampleVals(j,k)];
            xSampleIDsVec = [xSampleIDsVec,xSampleIDs{j,k}];
        end
    end
    plot(xPlot,yPlot,'o');
    xticks(xSampleValsVec);
    xticklabels(xSampleIDsVec);
    xtickangle(45);
    ylabel(geomMetricList{i});
end
