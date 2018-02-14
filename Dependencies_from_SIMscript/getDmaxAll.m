function DmaxAll=getDmaxAll(cellInfo)
% function to extract all the maximal Diameter from cellInfo and put
% them in a matrix for analysis.
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
DmaxAll=NaN(num_cells,num_frames);
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        try
            DmaxAll(cellIdx,frIdx)=cellInfo{frIdx,1}{1,cellIdx}.maxD;
        catch
        end
    end
end
end