function plGeom = loadPlasmGeom(folderName,nNuclei)
% load the Volume (in Voxels and microns) for each nucleoplasm and 
% stores is as a table where each row is a nucleus.
plGeom = [];
channelIdx =1;
for nucIdx = 1:nNuclei
    plasmIntFileName = fullfile(folderName,...
        ['C',num2str(channelIdx),'_nuc',num2str(nucIdx),'_plasmInt.csv']); 
    tmp = readtable(plasmIntFileName);
    varsToKeep = {'NumberOfVoxels','Volume'};
    tmp2 = tmp;
    for i=1:numel(tmp.Properties.VariableNames)
        if ~ismember(tmp.Properties.VariableNames{i},varsToKeep)
            tmp2 = removevars(tmp2,tmp.Properties.VariableNames{i});
        end
    end

    plGeom = [plGeom;tmp2];

end