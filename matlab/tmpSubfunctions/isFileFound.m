function [isFound,fileName] = isFileFound(folderName,filePattern)
% looks withing a folder (folderName) for a file containing the string
% filePattern
% returns isFOund =1 or 0 depending on whether file was found.
% returns fileName: the full name of the file (just one file, the last one found in the list)
    fileName = [];
    curFList = dir(folderName);
    isFound = 0;
    for j=1:numel(curFList)
        if contains(curFList(j).name,filePattern)
            fileName = fullfile(folderName,curFList(j).name);
            isFound = 1;
        end
    end

end