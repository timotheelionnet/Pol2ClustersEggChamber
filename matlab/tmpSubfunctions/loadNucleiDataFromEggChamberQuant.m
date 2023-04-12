folderName = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut2';
% expectation for folder architecture is:
% -input folder (arbitrary name)
%   - condition 1 (arbitrary name)
%       - eggChamber 1 (arbitrary name)
%           - "eggChamberCSV" 
%               - "allNucGeom.csv"
%               - "C1_allNucInt.csv"
%               - "C2_allNucInt.csv"
%                - ... for all channels
%               - "C1_wholeImgInt.csv"
%               - "C2_wholeImgInt.csv"
%                - ... for all channels
%               - "C1_eggChamberInt.csv"
%               - "C2_eggChamberInt.csv"
%                - ... for all channels
%       - eggChamber 2
%           - "eggChamberCSV" 
%   - condition 2
%       - eggChamber 1 (arbitrary name)
%           - "eggChamberCSV" 
%       - eggChamber 2
%           - "eggChamberCSV" 
% etc


% build eggChamberDataFolder object to simplify file loading etc
ec = eggChamberDataFolder(folderName);
ec.collectConditionsAndSamples();

% turn off the warning the Matlab modifies headers that contain
% unacceptable characters (like '.') to decluter the command line.
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%% load all nuclei data into a table nr
for i=1:ec.nConditions
    for j=1:ec.nConditions

        curNr = ec.loadEggChamberData(i,j);

    end
end