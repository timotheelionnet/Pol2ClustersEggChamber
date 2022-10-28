function g = findAllNucleiGroupCombinations()
% assigns nuclei 1 to 15 to four groups

% Gr 1: 4 largest nuclei, closest from oocyte
% Gr 2: 6 next largest
% Gr 3: 4 next largets
% Gr 4: 1 smallest, furthest from oocyte

% ouputs matrix g where:
% columns 1:4 correspond to nuclei assigned to group 1
% columns 5:10 correspond to nuclei assigned to group 2
% columns 11:14 correspond to nuclei assigned to group 3
% column 15 corresponds to nucleus assigned to group 4

n4 = 15;
idxP3 = nchoosek(1:14,4);
n3 = size(idxP3,1);
idxP2 = nchoosek(1:10,6);
n2 = size(idxP2,1);
nTot = n2*n3*n4;
g = zeros(nTot,15);

% indices for group 4
g(:,15) = reshape(repmat((1:n4),n2*n3,1),nTot,1);

% get indices for group 3
for i=1:n4
    r3 = setdiff(1:n4,i);
    g( ((i-1)*n3*n2+1):(i*n3*n2), 11:14 ) = repmat(r3(idxP3),n2,1);
end

% get indices for group 2
for i=1:n4
    for j=1:n3
        r2 = setdiff(1:n4,[g((i-1)*n3*n2+j,15),g((i-1)*n3*n2+j,11:14)]);
        idx = (i-1)*n2*n3 + (0:(n2-1))*n3 + j;
        g(idx,5:10) = r2(idxP2);
    end
end

A = repmat(1:n4,size(g,1),1);
B = g(:,5:15);
mask = all(bsxfun(@ne,A,permute(B,[1 3 2])),3);
At = A.'; %//'
g(:,1:4) = reshape(At(mask.'),[],size(A,1)).';

end
