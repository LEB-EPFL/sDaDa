% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Aster Vanhecke 
% 
% Contact:
% e-mail:
% aster.vanhecke@epfl.ch
% 
% paper mail:
% EPFL SB IPHYS LEB 
% BSP 428 (Cubotron UNIL) 
% Rte de la Sorge 
% CH-1015 Lausanne
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function DmaxAll=getDmaxAll(cellInfo)
% function to extract all the maximal Diameter from cellInfo and put
% them in a matrix for analysis.

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