% remove all the cluster values outside of the volume range


% collect all the variables from the cluster table except those that were
% lifted from the nucleus table
varsToExclude = ismember(ec.clustT.Properties.VariableNames,{'clust_Label'});
prefixList = {'nuc_','sampleROI_','wholeImg_','eggChamber_','eggChamber_','plasm_','nucleoli_'};
for i=1:numel(prefixList)
    varsToExclude = varsToExclude  | ...
        cellfun(@startsWith,...
        ec.clustT.Properties.VariableNames,repmat(prefixList(i),...
        size(ec.clustT.Properties.VariableNames)));
end

varsToExclude(ismember(ec.clustT.Properties.VariableNames,{'nuc_Label'})) = 0;
varsToExclude(ismember(ec.clustT.Properties.VariableNames,{'eggChamber_Idx'})) = 0;
varsToExclude(ismember(ec.clustT.Properties.VariableNames,{'eggChamber_Stage'})) = 0;

n3 = ec.clustT(:,~varsToExclude);

% average each metric across the nucleus
n3 = grpstats(n3,...
                ["cond_Idx","sample_Idx","sample_InputFileName","eggChamber_Idx","eggChamber_Stage","nuc_Label"],...
                ["mean","std"]);

% rename the number of member in each nucleus into nClusters
n3 = renamevars(n3,{'GroupCount'},{'nClusters'});

% remove the new row names created by the grpstats function
n3.Properties.RowNames = {}; 

% replace prefixes mean_ and std_ by nucAvgClustMinVol and nucStdClustMinVol

% join the averaged metrics with the nucleus table
n4 = join(ec.nucFullT,n3,"Keys",["cond_Idx","sample_Idx","sample_InputFileName","eggChamber_Idx","eggChamber_Stage","nuc_Label"]);
