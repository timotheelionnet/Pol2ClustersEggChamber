function [volumeMetric,meanGroupVol,stdGroupVol] = computeNucleiGroupVolumeMetrics(n)
% n is a nuclear table that should hold 15 rows, one for each nuclei (one egg chamber)
% with columns that include Centroid_X, Centroid_Y, Centroid_Z

% assumes that the indices of the nuclei assigned to each nurse cell
% are ranked by group:
% First four nuclei in group 1, next 6 group 2, next 4 group 3 and
% last nuclei in group 4.

% expected goup order
gIdx = {[1,2,3,4];...
    [5,6,7,8,9,10];...
    [11,12,13,14];...
    15};

meanGroupVol = zeros(1,4);
stdGroupVol = zeros(1,4);
for i=1:4
    meanGroupVol(i) = mean(n.Volume(gIdx{i}) );
    stdGroupVol(i) = std(n.Volume(gIdx{i}) );
end

volumeMetric = (meanGroupVol(1) - meanGroupVol(2))/(stdGroupVol(1) + stdGroupVol(2))...
    + (meanGroupVol(2) - meanGroupVol(3)) /(stdGroupVol(2) + stdGroupVol(3))...
    + (meanGroupVol(3)-meanGroupVol(4))/(stdGroupVol(3) + stdGroupVol(4));

end