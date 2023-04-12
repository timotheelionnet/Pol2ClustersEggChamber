function cl = addNucleiMetricsToClusterTable(nr,cl,eggChamberID,nucID,nc)
    
    varsToAdd = {...
        'nucVoxelCount',...
        'nucVolume',...
        };

    channelVarsToAdd = {...
        'Mean_eggChamberCorr',...
        'Median_eggChamberCorr'};

    % find row indices in cluster table that belong to current nucleus 
    idxC = (cl.eggChamberID == eggChamberID) & (cl.nucID == nucID);

    % compute number of clusters
    nClusters = sum( double( idxC ) );
    
    % find index of corresponding nucleus in nucleus table
    idxN = (nr.eggChamberID == eggChamberID) & (nr.nucID == nucID);
    
    newVarsToAdd = varsToAdd;
    for i=1:nc
        for j=1:numel(channelVarsToAdd)
            varsToAdd = [varsToAdd, ['nucC',num2str(i),'_',channelVarsToAdd{j}]];   
            newVarsToAdd = [newVarsToAdd, ['nucC',num2str(i),'_',channelVarsToAdd{j}]];   
        end
                
    end
    
    % compile all metrics to insert into one table
    for i=1:numel(varsToAdd)
        if i==1
            t = table(nr.(varsToAdd{1})(idxN),'variableNames',newVarsToAdd(1));
        else
            x = nr.(varsToAdd{i})(idxN);
            t = addvars(t,x,'newVariableNames',newVarsToAdd{i});
        end 
    end

    if ~ismember('nucVoxelCount',cl.Properties.VariableNames)
       % create an empty (zeros) table to insert all the
       % cluster-related variables into the nucleus table
       t = repmat(t,size(cl,1),1);
       idx = ~idxC;
       t{idx,:} = 0;
       for v=1:numel(t.Properties.VariableNames)
           cl = addvars(cl,t.(t.Properties.VariableNames{v}),...
               'NewVariableNames',t.Properties.VariableNames(v),...
               'Before','clustC1_Mean_plasmCorr');
       end
    else
        for v=1:numel(t.Properties.VariableNames)
           cl.(t.Properties.VariableNames{v})(idxC) = ...
               repmat(t.(t.Properties.VariableNames{v})(1),nClusters,1);
       end
    end

end