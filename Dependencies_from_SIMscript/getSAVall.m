function [SAVall, SAall, Vall]=getSAVall(cellInfo)
% function to extract all the surface area, volume and SA/V data from 
% cellInfo and put them in a cell array for analysis.
% Does rescaling from pxl to micrometer, pixSize=0.03 um
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
SAall=NaN(num_cells,num_frames);
Vall=NaN(num_cells,num_frames);
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        try
            SAall(cellIdx,frIdx) = cellInfo{frIdx,1}{1,cellIdx}.fwhmSurf .*0.03^2;
            Vall(cellIdx,frIdx)  = cellInfo{frIdx,1}{1,cellIdx}.fwhmVol .*0.03^3;
        catch
            warning(['Failed for cell ' num2str(cellIdx) ' in frame ' num2str(frIdx)])
        end
    end
end
SAVall=SAall./Vall;
end