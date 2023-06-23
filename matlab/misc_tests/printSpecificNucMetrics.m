cond = 1;
sample = 2;
nucIdx = [20,3,9,7];

b = zeros(numel(nucIdx),8);
for i=1:numel(nucIdx)
    imgIdx = find(  ec.nucFullT.cond_Idx == cond ...
                & ec.nucFullT.sample_Idx == sample ...
                & ec.nucFullT.nuc_Label == nucIdx(i));
 
    cIdx = 1; 
    a = [ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 2; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 3; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    cIdx = 4; 
    a = [a,ec.nucFullT.(['nucleoli_C',num2str(cIdx),'Mean_raw'])(imgIdx),ec.nucFullT.(['plasm_C',num2str(cIdx),'Mean_raw'])(imgIdx)];
    
    b(i,:) = a;

end

b