function t2 = assignNucleiToGroupsSingleEggChamber(t)

% if there are not 15 nuclei in the egg chamber, assign all nuclei to group 0
if size(t,1) ~= 15
    ecGroup = zeros(size(t,1),1);
    t2 = addvars(t,ecGroup,'NewVariableNames','ecGroup','After','nucLabel');
    return;
end

% define axis between centroid of 4 largest nuclei (c1) and furthest apart
% nucleus (c2)
[~,idxVol] = sortrows(t,'nucVolume','descend');
ecGroup = zeros(15,1);
for i=1:4
    ecGroup(idxVol(1:4),1) = 1;
end

xyz = [t.nucCentroid_X,t.nucCentroid_Y,t.nucCentroid_Z];
c1 = mean( xyz(idxVol(1:4),:) );

d = norm( repmat(c1,15,1) - xyz); 
d(idxVol(1:4)) = -Inf;
[~,idxMaxDist] = max(d);
ecGroup(idxMaxDist,1) = 4;

c2 = xyz(idxMaxDist,:);

% rank nuclei based on position along axis s
s = (xyz - repmat(c1,15,1))*(c2 - c1)';

s(idxVol(1:4)) = Inf;
s(idxMaxDist) = Inf;

[~,idxS] = sort(s,'ascend');
ecGroup(idxS(1:6)) = 2;
ecGroup(idxS(7:10)) = 3;

t2 = addvars(t,ecGroup,'NewVariableNames','ecGroup','After','nucLabel');


end