function cellInfo=changeTz(cellInfo,Tzs)
% Function to tell cellInfo when the FtsZ/FtsW is assembled/arrived. Input
% is cellInfo and Tzs, a vector containing the Tz (arrival/assembly FRAME)
% for each cell in order of cell number. The function changes the field
% 'Zassembled', which can be used by analyzeShapeDynamics for constraining Tc in
% the fits.
% Author: Aster Vanhecke

num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
if size(Tzs)== [1 num_cells]
elseif size(Tzs)== [num_cells 1]
else
    error('Tzs does not have the right size!')
end
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        if isfield(cellInfo{frIdx,1}{1,cellIdx},'Zassembled')
            cellInfo{frIdx,1}{1,cellIdx}.Zassembled= frIdx >= Tzs(cellIdx);
        end
    end
end