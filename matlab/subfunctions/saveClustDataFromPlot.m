%% save plot as table and figures
function saveClustDataFromPlot(fhFull,fhNuc,fullTable,avgNucClustTable,avgEcClustTable,avgCondClustTable, folderName,figName)
    saveas(fhFull,fullfile(folderName,[figName,'AllClust.fig']));
    try
        saveas(fhFull,fullfile(folderName,[figName,'AllClust.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',figName,'AllClust - likely permission issue.']);
    end

    saveas(fhNuc,fullfile(folderName,[figName,'AvgByNuc.fig']));
    try
        saveas(fhNuc,fullfile(folderName,[figName,'AvgByNuc.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',figName,'AvgByNuc - likely permission issue.']);
    end

    writetable(fullTable,fullfile(folderName,[figName,'Full.txt']),'Delimiter','\t');
    writetable(avgNucClustTable,fullfile(folderName,[figName,'AvgNuc.txt']),'Delimiter','\t');
    writetable(avgEcClustTable,fullfile(folderName,[figName,'AvgEc.txt']),'Delimiter','\t');
    writetable(avgCondClustTable,fullfile(folderName,[figName,'AvgCond.txt']),'Delimiter','\t');
end