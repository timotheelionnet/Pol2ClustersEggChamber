function [ctr,spr] = computeNucleiGroupCoordinates(n)
% Computes the centroids of each eggchamber nuclei subgroup
% n is a nuclear table that should hold 15 rows, one for each nuclei (one egg chamber)
% with columns that include Centroid_X, Centroid_Y, Centroid_Z

% OUTPUT
% ctr: centroid coordinates of all groups, formatted as a 3 x 4 matrix, where
% col 1 is [ gr1_centroid_X ; gr1_centroid_Y ; gr1_centroid_Z ]
% the second column same for group 2 etc.

% spr is the sqrt of the mean square distance, formatted as a 3 x 4 matrix, where
% col 1 is [ gr1_std_X ; gr1_std_Y ; gr1_std_Z ]
% the second column same for group 2 etc.

%%

% gIdx{i} contains the indices of the nuclei assigned to each nurse cell
% group. First four nuclei in group 1, next 6 group 2, next 4 group 3 and
% last nuclei in group 4.
gIdx = {[1,2,3,4];...
    [5,6,7,8,9,10];...
    [11,12,13,14];...
    15};

ctr = zeros(3,4);
spr = zeros(3,4);
for i=1:4
    ctr(1,i) = mean(n.Centroid_X(gIdx{i}));
    ctr(2,i) = mean(n.Centroid_Y(gIdx{i}));
    ctr(3,i) = mean(n.Centroid_Z(gIdx{i}));

    spr(1,i) = std(n.Centroid_X(gIdx{i}),1);
    spr(2,i) = std(n.Centroid_Y(gIdx{i}),1);
    spr(3,i) = std(n.Centroid_Z(gIdx{i}),1);

end