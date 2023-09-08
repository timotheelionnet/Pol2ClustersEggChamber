
%load random cluster values to double check table is correct
n = 10;
idx = floor(size(ec.clustT,1)*rand(n,1));
for i=1:n
    
    
    disp(['Idx: ',num2str(idx(i))]);
    curClust = ec.clustT(idx(i),1:8)

    % find cluster folder
    clustFolder = [ec.inFolder,filesep,ec.conditionNames{ec.clustT.cond_Idx(idx(i))},filesep,ec.sampleNames{ec.clustT.cond_Idx(idx(i))}{ec.clustT.sample_Idx(idx(i))},filesep,'nucCSV',filesep]

    %load nucleus geometry
    nucName = ['nuc',num2str(ec.clustT.nuc_Label(idx(i)))];
    r = readtable(fullfile(clustFolder,[nucName,'_localGeom.csv']));

    % compare nucleus volume across tables
    disp(['Nuc Volume from localGeom file: ',num2str(r.Volume)]);
    idxNuc = ec.nucFullT.cond_Idx == ec.clustT.cond_Idx(idx(i)) ...
        & ec.nucFullT.sample_Idx == ec.clustT.sample_Idx(idx(i)) ...
        & ec.nucFullT.nuc_Label == ec.clustT.nuc_Label(idx(i));
    v = ec.nucFullT.nuc_Volume( idxNuc );
    disp(['Nuc Volume from nuc table: ',num2str(v)]);
    
    v = ec.clustT.nuc_Volume( idx(i) );
    disp(['Nuc Volume from cluster table: ',num2str(v)]);
    
    disp(' ');

end