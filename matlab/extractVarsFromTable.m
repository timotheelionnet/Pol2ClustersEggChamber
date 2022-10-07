function nrs = extractVarsFromTable(nr,nv)

varsToRemove = setdiff(nr.Properties.VariableNames,nv);
        
nrs = nr;
for i=1:numel(varsToRemove)
    nrs = removevars(nrs,varsToRemove{i});

end