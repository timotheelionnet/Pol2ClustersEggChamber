function t2 = assignNucleiToGroupsWholeTable(t)

conds = unique(t.condIdx);
samples = unique(t.sampleIdx);

t2 = [];
% loop through eggchambers
for i=1:numel(conds)
    for j=1:numel(samples)
        curT = t(t.sampleIdx == samples(j) & t.condIdx == conds(i),:);
        t2 = [t2;assignNucleiToGroupsSingleEggChamber(curT)];
    end
end


end

