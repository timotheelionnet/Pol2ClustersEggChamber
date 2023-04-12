function [nc,allFilesPresent] = computeNumberOfChannelsInSubfolder(folderName)
% looks within a folder called folderName for files following the pattern
% C1-wholeNucleiMeas.csv, C2-wholeNucleiMeas.csv, etc

% based on the number of files found, returns nc the number of channels
% also double checks that all indices 1:nc are present and adjusts the
% value of the allFIlesPresent accordingly to 1 or 0.

% make sure expected file exists
fList = {};
curFList = dir(folderName);
for j=1:numel(curFList)
    if contains(curFList(j).name,'nucInt.csv')
        fList = [fList,curFList(j).name];
    end
end

nc = numel(fList);
allFilesPresent = 1;
for i=1:nc
    curName = ['C',num2str(i),'_nucInt.csv'];
    if ~ismember(curName,{curFList.name})
        disp(['While counting nuclei, could not find file ',curName,' in folder ',folderName]);
        allFilesPresent = 0;
    end
end

end