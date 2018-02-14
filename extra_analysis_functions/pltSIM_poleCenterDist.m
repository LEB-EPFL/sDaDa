% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Anna Archetti
% 
% Contact:
% e-mail:
% anna.archetti@epfl.ch
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

%------------------------------------------------------------------------%
% pltSIM_poleCenterDist
%------------------------------------------------------------------------%

function [p1C, p2C, t] = pltSIM_poleCenterDist(cellInfo, tInt, pixSize, divededVarEachFrm, tc)

    [frNum id]= size(cellInfo);
    [id cellNum] = size(cellInfo{1,1});
    
    for cellIdx = 1 : cellNum       
        for frIdx = 1 :  frNum      
            
            divededVarIdx = find(divededVarEachFrm(:,1) == frIdx);
            divededVarTemp = divededVarEachFrm(divededVarIdx, 2:3);
            if ~isempty(cellInfo{frIdx, 1}{1, cellIdx}) || divededVarTemp(cellIdx, 2) == 1
                
              
                if divededVarTemp(cellIdx, 2) == 0
                    c0 = cellInfo{frIdx, 1}{1, cellIdx}.c0;
                    leng = cellInfo{frIdx, 1}{1, cellIdx}.bactLength*pixSize; 
                    mesh = cellInfo{frIdx, 1}{1, cellIdx}.mesh; 
                    t = frIdx*tInt ;

                    
                    %       cell1  cell2 ... celln
                    % fr1   t1 w1  t2 w2 ... tn wn
                    % fr2
                    twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t c0];
                    tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t leng];
                else % if it is divided
                    
                     t = frIdx*tInt ;

                    %       cell1  cell2 ... celln
                    % fr1   t1 w1  t2 w2 ... tn wn
                    % fr2
                    twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t 0];
                    tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t nan];
                end
            else
                twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [nan nan];
                tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [nan nan];
            end
            
        end % end frame cycle
    end% end cell cycle
    

end