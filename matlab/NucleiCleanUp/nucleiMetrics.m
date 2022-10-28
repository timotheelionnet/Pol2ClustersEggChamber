%% filepaths
csvFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut";


%% load results into one big table

% collect subfolder names
sd = dir(csvFolder);
ctr = 0;
subDirList = {};
for i=1:numel(sd)
    if isfolder(fullfile(csvFolder,sd(i).name)) ... 
            && ~strcmp(sd(i).name,'.') && ~strcmp(sd(i).name,'..')
        ctr = ctr+1;
        subDirList{ctr,1} = fullfile(csvFolder,sd(i).name);
    end
end

% turn off the warning the Matlab modifies headers that contain
% unacceptable characters (like '.') to decluter the command line.
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% collect files
t = [];
for i=1:numel(subDirList)
    fl = dir(subDirList{i,1});
    for j=1:numel(fl)
        [~,~,e] = fileparts(fl(j).name);
        if strcmp(e,'.csv')
            curT = readtable(fullfile(subDirList{i,1},fl(j).name));
            % add conditionID column i.e. index corresponding to the
            % subfolder it came from 
            cond = i*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','conditionID'); 

            % add eggChamberID column i.e. index corresponding to the
            % file it came from 
            cond = j*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','eggChamberID');
            
            % append to global table
            t = [t;curT];
        end
    end
end

%% load table w annotations
t2 = readtable('/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut/MatlabTableFirstPassWAnnotation.txt');


%% perform PCA on table w/o annotations
t3 = removevars(t2,{'isWrong','Centroid_X','Centroid_Y','Centroid_Z','conditionID','Label','Var1'});
t4 = t3(:,vartype("numeric"));
t4 = table2array(removevars(t4,{'eggChamberID'}));
ids = unique(t3.eggChamberID);
t5 = [];
for i=1:numel(ids)
    idx = t3.eggChamberID == ids(i);
    t5 = [ t5; zscore(t4(idx,:))];
end

[coeff,score,latent] = pca(t5);

figure('Name','PCA'); hold;
plot(score(t2.isWrong==0,1),score(t2.isWrong==0,2),'o');
plot(score(t2.isWrong==1,1),score(t2.isWrong==1,2),'o');

%% plot stuff
figure; hold;
plot(t2.Volume(t2.isWrong==0),t2.MeanBreadth(t2.isWrong==0),'o');
plot(t2.Volume(t2.isWrong==1),t2.MeanBreadth(t2.isWrong==1),'o');

%%
figure; hold;
plot(t2.Volume(t2.isWrong==0),t2.dist(t2.isWrong==0),'o');
plot(t2.Volume(t2.isWrong==1),t2.dist(t2.isWrong==1),'o');

%%
figure; hold;
plot(t2.Volume(t2.isWrong==0),t2.SurfaceArea(t2.isWrong==0)./t2.Volume(t2.isWrong==0),'o');
plot(t2.Volume(t2.isWrong==1),t2.SurfaceArea(t2.isWrong==1)./t2.Volume(t2.isWrong==1),'o');



