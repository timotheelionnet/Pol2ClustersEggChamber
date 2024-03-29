%% use at your own risk!!! this script will delete tiff files within a root folder
rootFolder = '/Users/lionnt01/Dropbox/data/feiyue/20230508_PROTAC';
% rootFolder = '/Users/lionnt01/Dropbox/data/feiyue/nucSeg20_3img/out_v4';

% this flag when set to 1 will retain the 2D masks outlining the egg
% chambers. When set to zero, the 2D masks will be erased. 
% recommend 1
keepEggChamberMasks = 1;

% this flag when set to 1 will retain the final eggChamber tifs. when set
% to zero, those files will be erased.
% recommend 1
keepFinalEggChamberTif = 1;

% this flag when set to 1 will retain the init eggChamber tifs. when set
% to zero, those files will be erased.
% recommend 0
keepInitEggChamberTif = 0;

% this flag when set to 1 will retain the plasm-corrected nuclei tifs. when set
% to zero, those files will be erased.
% recommend 0
keepPlasmCorrNucTif = 0;

% this flag when set to 1 will retain the raw nuclei tifs. when set
% to zero, those files will be erased.
% recommend 1
keepRawNucTif = 1;

%% Collect all tifs in root folder
fileList = dir(fullfile(rootFolder, '**',filesep,'*.*'));  %get list of files and folders in any subfolder
fileList = fileList(~[fileList.isdir]);  %remove folders from list

% keep only tiff files
idx1 = cell2mat( ...
    cellfun( @contains, {fileList(:).name}, ...
    repmat({'tif'},size({fileList(:).name})),...
    repmat({'Ignorecase'},size({fileList(:).name})),...
    repmat({true},size({fileList(:).name})),...
    'UniformOutput',0 ))';  
idx2 = cell2mat( ...
    cellfun( @contains, {fileList(:).name}, ...
    repmat({'tiff'},size({fileList(:).name})),...
    repmat({'Ignorecase'},size({fileList(:).name})),...
    repmat({true},size({fileList(:).name})),...
    'UniformOutput',0 ))';  
idx1 = idx1 | idx2;

tiffFileList = fileList(idx1);
clear idx1 idx2 

%% colect list of files to delete
if ~keepEggChamberMasks
    idxToDeleteEggChamberMasks = ~cellfun(@isempty,cellfun( @regexp, {tiffFileList(:).name}, ...
        repmat({'_eggChamber\d*[.]tif$'},size({tiffFileList(:).name})),...
        'UniformOutput',0 ))';  
else
    idxToDeleteEggChamberMasks = false(size({tiffFileList(:).name}))';
end

if ~keepFinalEggChamberTif
    idxToDeleteFinalEggChamberTif = ~cellfun(@isempty,cellfun( @regexp, {tiffFileList(:).name}, ...
        repmat({'FinalNucMask[.]tif$'},size({tiffFileList(:).name})),...
        'UniformOutput',0 ))';  
else
    idxToDeleteFinalEggChamberTif = false(size({tiffFileList(:).name}))';
end

if ~keepInitEggChamberTif
    idxToDeleteInitEggChamberTif = ~cellfun(@isempty,cellfun( @regexp, {tiffFileList(:).name}, ...
        repmat({'InitNucMask[.]tif$'},size({tiffFileList(:).name})),...
        'UniformOutput',0 ))';  
else
    idxToDeleteInitEggChamberTif = false(size({tiffFileList(:).name}))';
end

if ~keepPlasmCorrNucTif
    idxToDeletePlasmCorrNucTif = ~cellfun(@isempty,cellfun( @regexp, {tiffFileList(:).name}, ...
        repmat({'nuc\d*_masks_plasmCorr[.]tif$'},size({tiffFileList(:).name})),...
        'UniformOutput',0 ))';  
else
    idxToDeletePlasmCorrNucTif = false(size({tiffFileList(:).name}))';
end

if ~keepRawNucTif
    idxToDeleteRawNucTif = ~cellfun(@isempty,cellfun( @regexp, {tiffFileList(:).name}, ...
        repmat({'nuc\d*_masks[.]tif$'},size({tiffFileList(:).name})),...
        'UniformOutput',0 ))';  
else
    idxToDeleteRawNucTif = false(size({tiffFileList(:).name}))';
end

disp(' ');
disp('******* LIST OF FILES TO DELETE *********');
filesToRemove = tiffFileList(   idxToDeleteEggChamberMasks ...
                                | idxToDeleteFinalEggChamberTif...
                                | idxToDeleteInitEggChamberTif...
                                | idxToDeletePlasmCorrNucTif...
                                | idxToDeleteRawNucTif);

for i=1:numel(filesToRemove)
    disp(fullfile(filesToRemove(i).folder,filesToRemove(i).name));
end

%% Delete Files
prompt = ['Do you wish to delete the '...
                ,num2str(numel(filesToRemove)),' files above? This is final Y/N [Y]: '];
txt = input(prompt,"s");
if isempty(txt)
    txt = 'N';
end
if strcmpi(txt,'Y')
    disp('Deleting files...');
    for i=1:numel(filesToRemove)
        delete(fullfile(filesToRemove(i).folder,filesToRemove(i).name));
    end
    disp('Done.');
else
    disp('No files deleted.');
end


