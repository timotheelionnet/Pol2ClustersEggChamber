function t2 = addAnnotationsAndVarsToNucleiTable(t,pathToAnnotation)

t2 = readtable(pathToAnnotation);

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

end