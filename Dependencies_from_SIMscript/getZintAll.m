function ZintAll=getZintAll(cellInfo)
% function to extract all the Z intensity profiles from cellInfo and put
% them in a cell array for analysis. E.g. kymographs
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
ZintAll{num_cells,num_frames}=[];
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        try
            ZintAll{cellIdx,frIdx}=cellInfo{frIdx,1}{1,cellIdx}.ZintLength;
        catch
        end
    end
end
end