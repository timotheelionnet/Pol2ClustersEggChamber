function pl = loadPlasmInt(folderName,channelIdx, nNuclei)

pl = [];
for nucIdx = 1:nNuclei
    plasmIntFileName = fullfile(folderName,...
        ['C',num2str(channelIdx),'_nuc',num2str(nucIdx),'_plasmInt.csv']); 
    tmp = readtable(plasmIntFileName);
    oldVars = {'Mean','StdDev','Max','Min',...
        'Median','Mode','Skewness','Kurtosis'};
    tmp2 = tmp;
    for i=1:numel(tmp.Properties.VariableNames)
        if ~ismember(tmp.Properties.VariableNames{i},oldVars)
            tmp2 = removevars(tmp2,tmp.Properties.VariableNames{i});
        end
    end

    pl = [pl;tmp2];
end

