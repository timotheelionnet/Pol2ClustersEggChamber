function [nChannels,c] = findChannelNumberFromNucTable(nrs)
 x = regexp(nrs.Properties.VariableNames,'nucC(\d*)_Mean_raw','tokens'); 
 x = x(~cellfun(@isempty,x)); 
 x = [x{:}]; 
 x = [x{:}]; 
 nChannels = numel(unique(x));
 c = cellfun(@str2double,unique(x));
 c = sort(c);