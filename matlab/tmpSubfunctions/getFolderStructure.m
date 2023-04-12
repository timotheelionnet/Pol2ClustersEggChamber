classdef eggChamberDataFolder < handle
% from the below expected folder architecture, collects conditions and
% samples
% fs.inFolder = inputFolder
% fs.conditions(1).name = 'condition1'
% fs.conditions(1).sampleList = {'eggChamber1','eggChamber2',...}
% -input folder (arbitrary name)
%   - condition1 (arbitrary name)
%       - eggChamber1 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%       - eggChamber2 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%   - condition2
%       - eggChamber1 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 
%       - eggChamber2 (arbitrary name)
%           - "eggChamberCSV" 
%           - "eggChamberTIF" 


properties (GetAccess = 'public', SetAccess = 'public')
end


% collect list of subfolders
d = dir(folderName);
dFolders = d([d(:).isdir]);
dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));

end