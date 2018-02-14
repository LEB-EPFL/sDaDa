function [cellInfo,divededVarEachFrm] = tooLateCorrect(cellInfo, divededVarEachFrm)
% Function to correct when cells where uncorrectly classified as not
% divided during semi-auto analysis.
% Author: Aster Vanhecke

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


