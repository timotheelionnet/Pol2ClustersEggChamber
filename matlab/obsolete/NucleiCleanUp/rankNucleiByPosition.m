function [vidx,maxDist,maxIdx] =  rankNucleiByPosition(n)
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

% find all possible choices of 4 out of 15 for the group 1 (largest nuclei)
idx = nchoosek(1:15,4);

% keep only the conditions where the min nuclear volume in group 1 is
% larger than the median group volume in the rest of the nuclei
vidx = [];
for i=1:size(idx,1)
    i1 = n.Label(idx(i,:)); % indices of nuclei in group 1
    i234 = setdiff(n.Label,i1); % indices of all other nuclei
    
    if min(n.Volume(i1)) > median(n.Volume(i234))
        vidx = [vidx; reshape(idx(i,:),1,4)];
    end
end

% for each plausible condition based on volume, compute the axial distance between the
% centroid of group 1 and each other nucleus. Keep max distance. 
% use the combination that maximizes the max distance. 
maxDist = zeros(size(vidx,1),1);
maxIdx = zeros(size(vidx,1),1);
for i=1:size(vidx,1)
    % compute centroid of group 1 nuclei
    i1 = n.Label(vidx(i,:));
    i234 = setdiff(n.Label,i1);

    c1 = mean( [n.Centroid_X(i1),n.Centroid_Y(i1),n.Centroid_Z(i1)]);
    
    curDist = zeros(numel(i234),1);
    for j=1:numel(i234)
        jj = n.Label == i234(j);
        curDist(j) = norm(...
            [n.Centroid_X(jj),n.Centroid_Y(jj),n.Centroid_Z(jj)] -c1 );
    end
    [maxDist(i),idxj] = max(curDist);
    maxIdx(i) = i234(idxj);
end


% order nuclei into groups based on their rank along the axial coordinate. 

