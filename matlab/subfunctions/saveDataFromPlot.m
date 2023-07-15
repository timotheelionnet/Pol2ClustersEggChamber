%% save plot as table and figures
function saveDataFromPlot(fh,fullTable,avgEcTable,avgCondTable, folderName,figName)
    saveas(fh,fullfile(folderName,[figName,'.fig']));
    saveas(fh,fullfile(folderName,[figName,'.eps']),'epsc');
    writetable(fullTable,fullfile(folderName,[figName,'Full.txt']),'Delimiter','\t');
    writetable(avgEcTable,fullfile(folderName,[figName,'AvgEc.txt']),'Delimiter','\t');
    writetable(avgCondTable,fullfile(folderName,[figName,'AvgCond.txt']),'Delimiter','\t');
end