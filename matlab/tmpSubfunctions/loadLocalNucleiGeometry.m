
function ng = loadLocalNucleiGeometry(ng,nNuclei,folderName)
% load local coordinates, i.e. relative to the cropped nuclei
% images. These will be used to transform cluster coordinates from
% the local reference (cropped image) to a global reference (egg
% chamber).

for i=1:nNuclei
    [geomFileFound,gf] = isFileFound(folderName,['nuc',num2str(i),'_localGeom.csv']);
    if geomFileFound == 0
        disp(['Could not find nucleus ',num2str(i),' local geometry file in subfolder ',folderName,'; excluding from analysis.']);
        ng = [];
        return
    else
        disp(['Loading nucleus ',num2str(i),' local geometry from ',folderName,'...']);
    end
    lg = readtable(gf);
    
    % replace Fiji variable names X, Y, Z with localX, localY, localZ
    % respectively
    oldVars = {'Box.X.Min','Box.X.Max','Box.Y.Min',...
        'Box.Y.Max','Box.Z.Min','Box.Z.Max',...
        'Centroid.X','Centroid.Y','Centroid.Z',...
        'Elli.Center.X','Elli.Center.Y','Elli.Center.Z'};
    oldVars = cellfun(@strrep,oldVars,...
        repmat({'.'},size(oldVars)),repmat({'_'},size(oldVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,oldVars,...
        repmat({'X'},size(oldVars)),repmat({'localX'},size(oldVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,newVars,...
        repmat({'Y'},size(newVars)),repmat({'localY'},size(newVars)), ...
        'UniformOutput',0);

    newVars = cellfun(@strrep,newVars,...
        repmat({'Z'},size(newVars)),repmat({'localZ'},size(newVars)), ...
        'UniformOutput',0);
    
    for j=1:numel(newVars)
        newVars{j} = ['nuc',newVars{j}];
        ng.(newVars{j})(i) = lg.(oldVars{j})(1);
    end

end


