function cellInfo=getSAVforCellInfo(cellInfo,smooth)
% Function to run surface area and volume calculations on all cells and
% frames in cellInfo. Just loops through cellInfo and runs getSurfAndVolFromCellInfo
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);

for cellIdx=1:num_cells
    for frIdx=1:num_frames
        try
            cellInfo{frIdx,1}{1,cellIdx}=getSurfAndVolFromCellInfo(cellInfo{frIdx,1}{1,cellIdx},smooth);
        catch
            warning(['Failed for cell ' num2str(cellIdx) ' in frame ' num2str(frIdx)])
        end
    end
end