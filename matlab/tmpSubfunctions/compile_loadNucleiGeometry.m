function ng = compile_loadNucleiGeometry(folderName)


geomFileFound = 0;
for j=1:numel(curFList
    if contains(curFList(j).name,'nucleiGeom.csv')
        gf = fullfile(folderName,curFList(j).name);
        geomFileFound = 1;
    end
end

if geomFileFound == 0
    disp(['Could not find nuclei geometry file in subfolder ',folderName,'; excluding from analysis.']);
    ng = [];
    return
else
    disp(['Loading nuclei geometry from ',folderName,'...']);
end

ng = readtable(gf);

end
