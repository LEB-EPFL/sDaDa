function diamFWHMAll=getDiamFWHMAll(cellInfo)
% function to extract all the diamFWHM profiles from cellInfo and put
% them in a cell array for analysis.
% E.g. kymographs: do rescale first:
% diamFWHMAllResc=cellfun(@(x) x.*0.03, diamFWHMAll,'UniformOutput',false);
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
diamFWHMAll{num_cells,num_frames}=[];
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        try
            diamFWHMAll{cellIdx,frIdx}=cellInfo{frIdx,1}{1,cellIdx}.diamFWHM;
        catch
        end
    end
end
end