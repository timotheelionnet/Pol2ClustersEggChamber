function rawInt = loadClusterRawInt(folderName,curFolder,channelID,nucleusID)

    fileName = ['C',num2str(channelID),'_nucleus',num2str(nucleusID),'-clusterResults.csv'];
    
    rawInt = readtable(fullfile(fullfile(folderName,curFolder),fileName));


end