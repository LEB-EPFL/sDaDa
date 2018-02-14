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