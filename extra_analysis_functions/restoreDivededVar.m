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

function divededVarEachFrm=restoreDivededVar(cellInfo)
 % function to restore divededVarEachFrm based on cellInfo in case you
 % accidentally changed or deleted it.
 % Warning: in some cases this function classifies cells as divided, i.e.
 % in frames where it was just deleted the frame(s) before division.
 
 num_frames=size(cellInfo,1);
 num_cells=size(cellInfo{1,1},2);
 divededVarEachFrm=zeros(num_cells*num_frames,3);
 for cellIdx=1:num_cells
     testDel=zeros(1,num_frames);
    for frIdx=1:num_frames
        divededVarEachFrm((frIdx-1)*num_cells+cellIdx,:)=[frIdx cellIdx isempty(cellInfo{frIdx,1}{1,cellIdx})];
        testDel(frIdx)=isempty(cellInfo{frIdx,1}{1,cellIdx});
    end
    % deleted(empty) does not mean divided!
    lastCellStanding=find(testDel==0,1,'last');
    if sum(testDel(1:lastCellStanding)) % if there was any cell deleted before "division"
        frIdx=find(testDel(1:lastCellStanding)==1);
        for frIdx=frIdx
            divededVarEachFrm((frIdx-1)*num_cells+cellIdx,:)=[frIdx cellIdx 0]; % these are not divided
        end
    end
 end
 