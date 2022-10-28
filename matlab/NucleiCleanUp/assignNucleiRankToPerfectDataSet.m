%% filepaths
csvFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut";
t = loadPreSegmentedNucGeomIntCSVs(csvFolder);

%% load table w annotations
% includes manual annotations that say which nuclei are correctly vs incorrectly segmented
% also computes extra variables, e.g. surfToVolRatio, DistRank (rank of the
% distance from center), CV (coeff of variation of the DAPI intensity)
pathToAnnotation = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut/MatlabTableFirstPassWAnnotation.txt';
t2 = addAnnotationsAndVarsToNucleiTable(t,pathToAnnotation);

%% keep only correctly segmented nuclei as per manual annotation
tc = t2(t2.isWrong==0,:);

%%
cds = unique(tc.conditionID);
nEC = 0;
for j=1:numel(cds)
    nEC = max(nEC, numel(unique(tc.eggChamberID(tc.conditionID==cds(j)))));
end
metrics = cell(numel(cds),nEC);
for j=1:numel(cds)
    ids = unique(tc.eggChamberID(tc.conditionID==cds(j)));
    for i=1:numel(ids)
        % display filename in command line window
        disp(' ');
        
        % collect n, table with current eggchamber nuclei only
        idx = tc.eggChamberID == ids(i) & tc.conditionID==cds(j);
        n = tc(idx,:);
        disp(n.inputFileName{1});

        % compute and plot metrics pertaining to nuclei size and relative
        % orientation
        [metrics{i,j},g] = ...
            plotNucleiPatterningMetricsForAllGroupCombinations(...
            n,n.inputFileName{1});

        if ~isempty(g)
            % display in command line top 5 of the nuclei orders
            [~,p] = sort(metrics{i,j}(:,end),'ascend');
            r = 1:length(metrics{i,j}(:,end));
            r(p) = r;
            r = r(1:5);
            for k =1:5
                curG = n.Label(g(r(k),:))';
                disp(['rank ',num2str(k),' order: ',num2str(curG(1:4)),...
                    ' ',num2str(curG(5:10)),' ',...
                    num2str(curG(11:14)),' ',num2str(curG(15))]);
            end
        end
    end
end