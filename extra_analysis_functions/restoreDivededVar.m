 function divededVarEachFrm=restoreDivededVar(cellInfo)
 % function to restore divededVarEachFrm based on cellInfo in case you
 % accidentally changed or deleted it.
 % Warning: in some cases this function classifies cells as divided, i.e.
 % in frames where it was just deleted the frame(s) before division.
 % Author: Aster Vanhecke
 
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
 