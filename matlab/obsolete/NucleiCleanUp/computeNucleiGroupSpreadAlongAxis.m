function [maxStdAxialSpr,meanStdAxialSpr,meanStdAxialSprNorm,totAxialSpr,deltaTMetric,meanAxialSpr,stdAxialSpr] = computeNucleiGroupSpreadAlongAxis(n,ptsl)
% n is a nuclear table that should hold 15 rows, one for each nuclei (one egg chamber)
% with columns that include Centroid_X, Centroid_Y, Centroid_Z

% assumes that the indices of the nuclei assigned to each nurse cell
% are ranked by group:
% First four nuclei in group 1, next 6 group 2, next 4 group 3 and
% last nuclei in group 4.

% OUTPUT

%%

% expected goup order
gIdx = {[1,2,3,4];...
    [5,6,7,8,9,10];...
    [11,12,13,14];...
    15};

%format nuclei centroid coordinates as a 3 x n array - each column is [x;y;z]
pts = [n.Centroid_X';n.Centroid_Y';n.Centroid_Z'];

% compute coordinates of the projection of each nucleus along the line axis
t = computeParameterOfLinProjectionAlongAxis(pts,ptsl);

% total spread along axis of nuclei
totAxialSpr = max(t) - min(t);

meanAxialSpr = zeros(1,4); % mean axis coordinate for each nuclei group
stdAxialSpr = zeros(1,4); % std of the axis coordinate for each nuclei group
msd = zeros(1,4);
for i=1:4
    meanAxialSpr(1,i) = mean(t(gIdx{i}));
    stdAxialSpr(1,i) = std(t(gIdx{i}),1);

    p = pts(:,gIdx{i});
    msd(1,i) =  sqrt(sum(norm(p - mean(p,2)).^2,2))/size(p,2);
end

% max spread of the axis coordinate within a group
maxStdAxialSpr = max(stdAxialSpr);

% mean spread of the axis coordinate within a group
meanStdAxialSpr = mean(stdAxialSpr);

% mean spread of the axis coordinate divided by mean square distance within
% a group, i.e. fraction of the MSD explained by axial coordinate - should
% be low if the nuclei of the group have a similar position longitudinally.
x = stdAxialSpr./msd;
x(4) = 0.5; % arbitrary value to avoid NaN
meanStdAxialSprNorm = mean(x);


%% compute metric based on deltaT
% where deltaT is the axis coordinate difference between adjacent groups.
% CVdeltaT should be very low if groups are sorted along the axis, with similar
% deltaT between each neighbours. If groups are randomly ordered, positive
% and negative deltaTs should cancel out, i.e. deltaTMetric should be very small.
deltaTMetric = (meanAxialSpr(2) - meanAxialSpr(1))/(stdAxialSpr(1) + stdAxialSpr(2)) ...
    + (meanAxialSpr(3) - meanAxialSpr(2))/(stdAxialSpr(2) + stdAxialSpr(3)) ...
    + (meanAxialSpr(4) - meanAxialSpr(3))/(stdAxialSpr(3) + stdAxialSpr(4));

% using the difference between the most extreme groups as the orienting
% direction on the eggChamber axis
deltaTMetric = sign(meanAxialSpr(4) - meanAxialSpr(1)) * deltaTMetric;
