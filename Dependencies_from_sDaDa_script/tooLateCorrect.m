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

function [cellInfo,divededVarEachFrm] = tooLateCorrect(cellInfo, divededVarEachFrm)
% Function to correct when cells where uncorrectly classified as not
% divided during semi-auto analysis.

% get input
inputGood=false;
while ~inputGood
    clear cellIdx frameDivided
    cellIdx=input('Give cell numbers for cells that were not properly assigned as divided (put in [] for multiple cells):\n');
    if isempty('cellidx')
        return
    end
    frameDivided=input('Give the resp. frame numbers where they actually where divided (put in [] for multiple cells):\n');
    
    if size(cellIdx)==size(frameDivided)
        inputGood=true;
    end
end

% For all the cells that need correcting
for ij=1:length(cellIdx);
    % Correct divededVarEachFrm:
    isRightCell=divededVarEachFrm(:,2)==cellIdx(ij);
    isRightFrame=divededVarEachFrm(:,1)>=frameDivided(ij);

    framesCellWasDivided=zeros(size(isRightCell));
    for ii=1:length(isRightCell)
        framesCellWasDivided(ii)=isRightCell(ii) && isRightFrame(ii);
    end
    divededVarEachFrm(isRightCell,3)=framesCellWasDivided(isRightCell);
    
    % Correct cellInfo
    for frIdx=frameDivided(ij):size(cellInfo,1)
    cellInfo{frIdx,1}{1,cellIdx(ij)}=[];
    end
end


