
inputFolder = '/Users/lionnt01/Dropbox/data/feiyue/20230508_PROTAC';

disp('Warning: This script will not fix directory names');


fList = dir(fullfile(inputFolder, ['**',filesep,'*.*']));
chng = {};
k=0;
for i=1:numel(fList)
    curName = fullfile(fList(i).folder,fList(i).name);

    if ~fList(i).isdir && contains(fList(i).name,' ')
        oldName = fList(i).name;
        newName = strrep(fList(i).name,' ','_');
        if ~exist(fullfile(fList(i).folder,newName),'file')
            disp(['renaming ',fullfile(fList(i).folder,oldName),...
                'into ',fullfile(fList(i).folder,newName)]);
            movefile(fullfile(fList(i).folder,oldName),fullfile(fList(i).folder,newName));
               k = k+1;
               chng{k,1} = fullfile(fList(i).folder,oldName);
               chng{k,2} = fullfile(fList(i).folder,newName);               
        else
            disp(['WARNING! Could not rename file ',fullfile(fList(i).folder,oldName),...
                '; new file name ',fullfile(fList(i).folder,newName),' already exists']);
        end
    end  
end