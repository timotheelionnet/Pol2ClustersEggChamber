function [metrics,g] =  plotNucleiPatterningMetricsForAllGroupCombinations(n,figName)
addpath('../cbrewer');
% making sure nuclei group has 15 members
if size(n,1) ~= 15
    disp(['n has ',num2str(size(n,1)),' nuclei, stopping.']);
    metrics = [];
    g = [];
    return
end

%rank by decreasing volume
n = sortrows(n,'Volume','descend');

%g = findAllNucleiGroupCombinations();
g = exploreGroupCombinationsAroundBaseline();

% for i=1:size(g,1)
nc = size(g,1);
err = zeros(nc,1);
elong= zeros(nc,1);
maxStdAxialSpr = zeros(nc,1);
meanStdAxialSpr= zeros(nc,1);
meanStdAxialSprNorm= zeros(nc,1);
totAxialSpr= zeros(nc,1);
CVdeltaT= zeros(nc,1);
volumeMetric = zeros(nc,1);
deltaTMetric = zeros(nc,1);

for i=1:nc
        
    % compute group center coordinates and spread
    curN = n(g(i,:),:);
    [ctr,spr] = computeNucleiGroupCoordinates(curN);

    % linear fit of the group centers
    [ptsl,err(i),ctrf,dist,elong(i)] = linearFit3D(ctr);
    
    % extract metrics of how well the nuclei groups align along the egg
    % Chamber axis
    [maxStdAxialSpr(i),meanStdAxialSpr(i),meanStdAxialSprNorm(i),totAxialSpr(i),...
        deltaTMetric(i),~,~] = ...
        computeNucleiGroupSpreadAlongAxis(curN,ptsl);
    
    [volumeMetric(i),~] = computeNucleiGroupVolumeMetrics(curN);
    % compute sqrt(mean square spread along egg chamber axis) for each
    % group
%     figure; hold;
%     plot3(ctr(1,:),ctr(2,:),ctr(3,:),'o');
%     plot3(ctrf(1,:),ctrf(2,:),ctrf(3,:),'-');
%     plot3(curN.Centroid_X(1:4),curN.Centroid_Y(1:4),curN.Centroid_Z(1:4),'x');
%     plot3(curN.Centroid_X(5:10),curN.Centroid_Y(5:10),curN.Centroid_Z(5:10),'x');
%     plot3(curN.Centroid_X(11:14),curN.Centroid_Y(11:14),curN.Centroid_Z(11:14),'x');
%     plot3(curN.Centroid_X(15),curN.Centroid_Y(15),curN.Centroid_Z(15),'x');
%     axis equal
end

metrics = [ elong,...
            totAxialSpr,...
            volumeMetric,...
            deltaTMetric,...
            err,...
            maxStdAxialSpr./totAxialSpr, ...
            meanStdAxialSpr./totAxialSpr,...
            meanStdAxialSprNorm, ...
            ];

sortDir = { 'descend',...
            'descend',...
            'descend',...
            'descend',...
            'ascend',...
            'ascend',...
            'ascend',...
            'ascend',...
            };

rank =zeros(size(metrics));
for i=1:size(metrics,2)
    [~,idx] = sort(metrics(:,i),sortDir{i});
    idx = sortrows([(1:size(metrics,1))',idx],2,'ascend');
    rank(:,i) = idx(:,1);
end

metrics = [metrics,mean(rank,2)];

metricsNorm = sortrows(zscore(metrics),size(metrics,2),'ascend') ;

figure('Name',figName); 
ah1 = axes('Position',[0.3,0.55,0.65,0.35]);
ah2 = axes('Position',[0.3,0.1,0.65,0.35]);

imagesc(ah1,metricsNorm');
cm = cbrewer('div','RdYlBu',255);
colormap(ah1,flipud(cm));
colormapLims = [-2,2];
clim(ah1,colormapLims);
metricNames = {'elongation',...
    'total axial Spread',...
    'Normalize Vol diff between grps',...
    'axial coord diff across neighbors',...
    'grp ctr dist to fit',...
    'max grp axialSpr / totSpr',...
    'mean grp axialSpr / totSpr',...
    'mean grp axialSpr rel to grp. MSD',...
    'rank'};

yticklabels(ah1,metricNames);
yticks(ah1,1:size(metrics,2));
colorbar(ah1);

imagesc(ah2,metricsNorm(1:30,:)');
colormap(ah2,flipud(cm));
clim(ah2,colormapLims);
yticklabels(ah2,metricNames);
yticks(ah2,1:size(metrics,2));
colorbar(ah2);
%plot(elong,err,'o');
