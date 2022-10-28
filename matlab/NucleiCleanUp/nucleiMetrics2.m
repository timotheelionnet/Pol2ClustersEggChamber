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
        [~,f,e] = fileparts(fl(j).name);
        if strcmp(e,'.csv') && contains(f,'Geom')
            curT = readtable(fullfile(subDirList{i,1},fl(j).name));
            % add conditionID column i.e. index corresponding to the
            % subfolder it came from 
            cond = i*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','conditionID'); 

            % add eggChamberID column i.e. index corresponding to the
            % file it came from 
            cond = j*ones(size(curT,1),1);
            curT = addvars(curT,cond,'NewVariableNames','eggChamberID');
            
            %open matching intensity file
            curT2 = readtable(fullfile(subDirList{i,1},strrep(fl(j).name,'Geom','Int')));
            
            curT2 = removevars(curT2,...
                intersect(curT.Properties.VariableNames, curT2.Properties.VariableNames));
            % append to global table
            t = [t;[curT,curT2]];
        end
    end
end

%% load table w annotations
t2 = readtable('/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut/MatlabTableFirstPassWAnnotation.txt');

t2 = addvars(t,t2.isWrong,'NewVariableNames','isWrong');

% add surf to vol ratio
t2 = addvars(t2,t2.SurfaceArea./t2.Volume,'NewVariableNames','SurfToVolRatio');

% add distance rank
cds = unique(t2.conditionID);
distRank = zeros(size(t2,1),1);
for j=1:numel(cds)
    ids = unique(t2.eggChamberID(t2.conditionID==cds(j)));
    for i=1:numel(ids)
        idx = t2.eggChamberID == ids(i) & t2.conditionID==cds(j);
        curDist = t2.dist(idx);
        [~,curRank] = sort(curDist);
        distRank(idx) = curRank;
    end
end
t2 = addvars(t2,distRank,'NewVariableNames','distRank');

% add coeff of variation of the intensity distribution
t2 = addvars(t2,t2.StdDev./t2.Mean,'NewVariableNames','CV');

%% prepare t3, streamlined version of t2 for PCA

% remove non-numeric variables
t3 = t2(:,vartype("numeric"));

% remove variables that depend on the image size (absolute coordinates)
keepCentroidsForPCA = 0;
if ~keepCentroidsForPCA
    t3 = removevars(t3,{'Centroid_X','Centroid_Y','Centroid_Z',...
        'Box_X_Min','Box_X_Max','Box_Y_Min','Box_Y_Max','Box_Z_Min','Box_Z_Max',...
        'Elli_Center_X','Elli_Center_Y','Elli_Center_Z',...
        'InscrBall_Center_X','InscrBall_Center_Y','InscrBall_Center_X'});
end

% remove ID variables - keeping the label column since it tends to be
% biased towards bad nuclei at the end
t3 = removevars(t3,{'Var1','conditionID','isWrong','eggChamberID'});

% remove variables which include NaNs
idx = false(1,size(t3,2));
for i=1:size(t3,2)
    idx(i) = sum( isnan(t3.(t3.Properties.VariableNames{i})(:)) ) == 0;
end
t3 = t3(:,idx);

t3pca = table2array(t3);

%pcaZScoreMode = 'byEggChamber';
pcaZScoreMode = 'global';

switch pcaZScoreMode
    case 'global'
        t3z = zscore(t3pca);
    case 'byEggChamber'
        cds = unique(t2.conditionID);
        t3z = zeros(size(t3pca));
        for j=1:numel(cds)
            ids = unique(t2.eggChamberID(t2.conditionID==cds(j)));
            for i=1:numel(ids)
                idx = t2.eggChamberID == ids(i) & t2.conditionID==cds(j);
                t3z(idx,:) = zscore(t3pca(idx,:));
            end
        end
end
[coeff,score,latent] = pca(t3z);
figure('Name','PCA'); hold;


scatter3(score(t2.isWrong==0,1),score(t2.isWrong==0,2),score(t2.isWrong==0,3),'o');
scatter3(score(t2.isWrong==1,1),score(t2.isWrong==1,2),score(t2.isWrong==1,3),'o');

%% clean up based on PCA
t2clean = t2(score(:,1)<7 & score(:,1)>-10 & score(:,2)<4,:);
cds = unique(t2clean.conditionID);
t2clean2 = [];
for j=1:numel(cds)
    ids = unique(t2clean.eggChamberID(t2clean.conditionID==cds(j)));
    for i=1:numel(ids)
        idx = t2clean.eggChamberID == ids(i) & t2clean.conditionID==cds(j);
        curT2clean = sortrows(t2clean(idx,:),'Label','ascend');
        nRows = min(size(curT2clean,1),15);
        curT2clean = curT2clean(1:nRows,:);
        t2clean2 = [t2clean2;curT2clean];
    end
end

n0 = sum(t2.isWrong == 0);
n1 = sum(t2.isWrong == 1);
disp(['initially ',num2str(n0),' good nuclei and ',num2str(n1),' bad ones.']);

n0 = sum(t2clean.isWrong == 0);
n1 = sum(t2clean.isWrong == 1);
disp(['After threshold based on PCA, ',num2str(n0),' good nuclei and ',num2str(n1),' bad ones.']);

n0 = sum(t2clean2.isWrong==0);
n1 = sum(t2clean2.isWrong==1);
disp(['After limiting to 15 nuclei per egg chamber, ',num2str(n0),' good nuclei and ',num2str(n1),' bad ones.']);
%%
v = t2.Kurtosis;
figure; hold; plot(v(t2.isWrong==0),'o'); plot(v(t2.isWrong==1),'o');
%% cleanup based on mutliple criteria
maxVolume = 5000;
%minVolume = 200; %100 for less stringent
minVolume = 100; %this is the best, 200 for more stringent
%distRank weak separator by itself
maxSurfToVolRatio = 1; %this is the best, 0.8 for more stringent
%maxSurfToVolRatio = 0.8; %1 for less stringent
minSphericity = 0.4;
%maxLabel = 15; %17 for less stringent
maxLabel = 17; %this is the best, 15 for more stringent
maxLabel = 50; %for higher recall, slightly lower precision
minMeanBreadth = 7;
maxMeanBreadth = 25;
maxSkewness = 1;
maxKurtosis = 3;
maxCV = 0.45;
minCV = 0.15;

idxClean = t2.Volume < maxVolume & t2.Volume > minVolume ...
    & t2.SurfToVolRatio < maxSurfToVolRatio ...
    & t2.Sphericity > minSphericity ...
    & t2.Label < maxLabel ...
    & t2.MeanBreadth > minMeanBreadth & t2.MeanBreadth < maxMeanBreadth ...
    & t2.Skewness < maxSkewness ...
    & t2.Kurtosis < maxKurtosis ...
    & t2.CV < maxCV & t2.CV > minCV;

t2clean = t2(idxClean,:);

n0 = sum(t2.isWrong == 0);
n1 = sum(t2.isWrong == 1);
disp(['initially ',num2str(n0),' good nuclei and ',num2str(n1),' bad ones.']);

n0 = sum(t2clean.isWrong == 0);
n1 = sum(t2clean.isWrong == 1);
disp(['After threshold based on individual metrics, ',...
    num2str(n0),' good nuclei and ',num2str(n1),' bad ones']);
disp(['       i.e. ',...
    num2str(n0/sum(t2.isWrong == 0)),' recall; ',...
    num2str(n0/(n0+n1)),' precision']);

