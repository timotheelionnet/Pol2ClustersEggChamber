function ng = loadGlobalNucleiGeometry(folderName)
% look within a folder for a file called nucleiGeom.csv
% if found, load the file as a table ng.

% make sure expected file exists
[geomFileFound,gf] = isFileFound(folderName,'nucGeom.csv');
folderName = char(folderName);
if geomFileFound == 0
    disp(['Could not find nuclei geometry file in subfolder ',folderName,'; excluding from analysis.']);
    ng = [];
    return
else
    disp(['Loading nuclei geometry from ',folderName,'...']);
end

ng = readtable(gf);

ng = renamevars(ng,'Label','ID');

oldVars = ng.Properties.VariableNames;
newVars = oldVars;
for i=1:numel(oldVars)
    ng = renamevars(ng,oldVars{i},['nuc',oldVars{i}]);
end
