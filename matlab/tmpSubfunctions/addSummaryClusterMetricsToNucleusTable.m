function nr = addSummaryClusterMetricsToNucleusTable(nr,cl,eggChamberIdx,nucIdx,nc,corrMode)
    
    eggChamberIDs = unique(nr.eggChamberID);
    nucIDs = unique(nr.nucID);

    % find row indices in cluster table that belong to current nucleus 
    idxC = (cl.eggChamberID == ...
        eggChamberIDs(eggChamberIdx)) & (cl.nucID == nucIDs(nucIdx));
    
    % compute number of clusters
    nClusters = sum( double( idxC ) );
    
    t = table(nClusters,'VariableNames',{'nClusters'});

    if nClusters == 0
        t.clustMeanVolume(1) = 0;
        t.clustStdDevVolume(1) = 0;
        t.clustMeanNumberOfVoxels(1) = 0;
        t.clustStdDevNumberOfVoxels(1) = 0;
        for k = 1:nc
            varName = ['clustC',num2str(k),'_MeanMean_plasmCorr'];
            t.(varName)(1) = 0;
            varName = ['clustC',num2str(k),'_StdDevMean_plasmCorr'];
            t.(varName)(1) = 0;
            
            plasmFileName = ['C',num2str(k),'_nuc',...
                num2str(nucIDs(nucIdx)),'_plasmInt.csv'];
            
            idxN = (nr.eggChamberID == ...
                        eggChamberIDs(eggChamberIdx)) & (nr.nucID == nucIDs(nucIdx));

            csvFolder = fullfile(...
                fullfile( nr.rootFolder(idxN),nr.eggChamberFolder(idxN) ),...
                'csvFIJI');
            pl = readtable(fullfile(csvFolder,plasmFileName) );

            if k==1
                t.plasmNumberOfVoxels(1) = pl.NumberOfVoxels(1);
                t.plasmVolume(1) = pl.Volume(1);
            end

            varName = ['plasmC',num2str(k),'_Mean'];
            t.(varName)(1) = pl.Mean(1);

            varName = ['plasmC',num2str(k),'_StdDev'];
            t.(varName)(1) = pl.StdDev(1);
            
            varName = ['plasmC',num2str(k),'_Median'];
            t.(varName)(1) = pl.Median(1);
            
            wholeImgVarName = ['wholeImgC',num2str(k),'_',corrMode];
            varName = ['plasmC',num2str(k),'_Mean_wholeImgCorr'];
            t.(varName)(1) = pl.Mean(1) - nr.(wholeImgVarName)(idxN);
            
            varName = ['plasmC',num2str(k),'_Median_wholeImgCorr'];
            t.(varName)(1) = pl.Median(1) - nr.(wholeImgVarName)(idxN);

            eggChamberVarName = ['eggChamberC',num2str(k),'_',corrMode];
            varName = ['plasmC',num2str(k),'_Mean_eggChamberCorr'];
            t.(varName)(1) = pl.Mean(1) - nr.(eggChamberVarName)(idxN);
            
            varName = ['plasmC',num2str(k),'_Median_eggChamberCorr'];
            t.(varName)(1) = pl.Median(1) - nr.(eggChamberVarName)(idxN);
        end
    else
    
        % compute average cluster volume per nucleus
        t.clustMeanVolume(1) = mean(cl.clustVolume(idxC));

        % compute stdDev cluster volume per nucleus
        t.clustStdDevVolume(1) = std(cl.clustVolume(idxC));

        % compute average cluster NumberOfVoxels per nucleus
        t.clustMeanNumberOfVoxels(1) = mean(cl.clustVoxelCount(idxC));

        % compute stdDev cluster NumberOfVoxels per nucleus
        t.clustStdDevNumberOfVoxels(1) = std(cl.clustVoxelCount(idxC));
        
        % find the row of the nucleus table that corresponds to the current nucleus
        idxN = (nr.eggChamberID == ...
            eggChamberIDs(eggChamberIdx)) & (nr.nucID == nucIDs(nucIdx));

        % loop through channels
        for k = 1:nc
            % compute average cluster intensity per nucleus for each channel
            varName1 = ['clustC',num2str(k),'_Mean_plasmCorr'];
            varName2 = ['clustC',num2str(k),'_MeanMean_plasmCorr'];
            t.(varName2)(1) = mean(cl.(varName1)(idxC));
    
            % compute stdDev cluster intensity per nucleus for each channel
            varName2 = ['clustC',num2str(k),'_StdDevMean_plasmCorr'];
            t.(varName2)(1) = std(cl.(varName1)(idxC));

        end
    end

    if ~ismember('clustMeanVolume',nr.Properties.VariableNames)
       % create an empty (zeros) table to insert all the
       % cluster-related variables into the nucleus table
       t = repmat(t,size(nr,1),1);
       idxC = ...
           (nr.eggChamberID == eggChamberIDs(eggChamberIdx)) & (nr.nucID == nucIDs(nucIdx));
       idxC = ~idxC;
       t{idxC,:}=0;
       for v=1:numel(t.Properties.VariableNames)
           nr = addvars(nr,t.(t.Properties.VariableNames{v}),...
               'NewVariableNames',t.Properties.VariableNames(v),...
               'Before','nucC1_Mean_eggChamberCorr');
       end
    else
        idxC = ...
            (nr.eggChamberID == eggChamberIDs(eggChamberIdx)) & (nr.nucID == nucIDs(nucIdx));
        for v=1:numel(t.Properties.VariableNames)
           nr.(t.Properties.VariableNames{v})(idxC) = ...
               t.(t.Properties.VariableNames{v})(1);
       end
    end
end