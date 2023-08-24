%% save plot as table and figures
function saveNucDataFromPlot(fh,fullTable,avgEcTable,avgCondTable, folderName,figName)
    saveas(fh,fullfile(folderName,[figName,'.fig']));
    try
        saveas(fh,fullfile(folderName,[figName,'.eps']),'epsc');
    catch 
        disp(['Could not save eps file for ',figName,' - likely permission issue.']);
    end
    writetable(fullTable,fullfile(folderName,[figName,'Full.txt']),'Delimiter','\t');
    writetable(avgEcTable,fullfile(folderName,[figName,'AvgEc.txt']),'Delimiter','\t');
    writetable(avgCondTable,fullfile(folderName,[figName,'AvgCond.txt']),'Delimiter','\t');
end