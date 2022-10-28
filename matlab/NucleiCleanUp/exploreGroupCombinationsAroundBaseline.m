function g = exploreGroupCombinationsAroundBaseline()
% assigns nuclei 1 to 15 to four groups using minor changes from
% baseline order (assumed to be 1 to 15)

% Gr 1: 4 largest nuclei, closest from oocyte
% Gr 2: 6 next largest
% Gr 3: 4 next largets
% Gr 4: 1 smallest, furthest from oocyte

% ouputs matrix g where:
% columns 1:4 correspond to nuclei assigned to group 1
% columns 5:10 correspond to nuclei assigned to group 2
% columns 11:14 correspond to nuclei assigned to group 3
% column 15 corresponds to nucleus assigned to group 4

g0 = 1:15;
n = [4,6,4,1];
r{1} = 1:4;
r{2} = 5:10;
r{3} = 11:14;
r{4} = 15;

% compute all single-pair permutations that switch nuclei between groups
g = g0;
for i=1:3
    for j=i+1:4
        % permutation matrix
        pm = [repmat((1:n(i))',n(j),1), ...
            reshape(repmat(1:n(j),n(i),1),n(i)*n(j),1)];
        
        % extract corresponding ranks before permutation
        r1 = reshape(r{i}(pm(:,1))',n(i)*n(j),1);
        r2 = reshape(r{j}(pm(:,2))',n(i)*n(j),1);
        rm = [r1,r2];

        %generate a matrix with the baseline order and replace ranks
        curG = repmat(g0,size(rm,1),1);
        for k=1:size(rm,1)
            curG(k,rm(k,1)) = g0(rm(k,2));
            curG(k,rm(k,2)) = g0(rm(k,1));
        end
        g = [g;curG];
    end
end

% add all permutations of the nuclei which maintain order
for iStart = 1:14
    g = [g;g0(mod(iStart+(1:15),15)+1)];
end

% add all permutation where one entry is removed, the order of the following
% nuclei are shifted left by one, and the removed entry is placed elsewhere
for i=1:15
    a = [1:i-1,i+1:15];
    if i~=1
        g = [g; [i,a]];
    end
    if i~=15
        g = [g; [a,i]];
    end
    for j=1:13
        g = [g; [ a(1:j),i,a(j+1:14) ]];
    end

end

% compute all two-pair permutations that switch nuclei between 2 groups
% ie. nuclei a and b of group 1 exchange w nuclei c and d of group 2

for i=1:3
    for j=(i+1):4
        if n(i)>1 && n(j)>1
            pi = [];
            for i1 = 1:n(i)
                for i2 = (i1+1):n(i) 
                    pi = [pi;i1,i2];
                end
            end
    
            pj = [];
            for j1 = 1:n(j)
                for j2 = (j1+1):n(j) 
                    pj = [pj;j1,j2];
                end
            end
            

            pm = [];
            for i1=1:size(pi,1)
                pm = [pm; ...
                    repmat(pi(i1,1:2),size(pj,1),1),pj];
            end
    
            % extract corresponding ranks before permutation
            r1 = r{i}(pm(:,1:2));
            r2 = r{j}(pm(:,3:4));
    
            % permutation matrix
            rm = [r1,r2];
    
            % generate a matrix with the baseline order and replace ranks
            curG = repmat(g0,size(rm,1),1);
            for k=1:size(rm,1)
                curG(k,rm(k,1)) = g0(rm(k,3));
                curG(k,rm(k,2)) = g0(rm(k,4));
                curG(k,rm(k,3)) = g0(rm(k,1));
                curG(k,rm(k,4)) = g0(rm(k,2));
            end
            g = [g;curG];
        end
    end
end


