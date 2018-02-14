function ElMat=getElMat(tlMatrix)
% get matrix with kL and L0 from fitting L=L0*exp(kL*T) to tlMatrix
num_cells=size(tlMatrix,2)/2;
ElMat(num_cells,2)=0;
for cellIdx=1:num_cells
    try
        [kL, l0]  = fitLengthVsTime(tlMatrix( ~isnan(tlMatrix(:,cellIdx*2)) ,cellIdx*2-1), tlMatrix(~isnan(tlMatrix(:,cellIdx*2)),cellIdx*2)); % fit t and l, but ignore nan's
    catch
        kL=NaN; l0=NaN;
    end
    ElMat(cellIdx,:)=[kL, l0];
end

% function ElMat=getElMatFromCellInfo(cellInfo)
% % function to create elongation matrix
% % takes L(t) for each cell in cell cellInfo, fits to an exponential and
% % puts the fit parameters (L(0) and exponent kL) in a matrix
% % Author: Aster Vanhecke
% num_frames=size(cellInfo,1);
% num_cells=size(cellInfo{1,1},2);
% ElMat=zeros(num_cells,2);
% for cellIdx=1:num_cells
%     L=zeros(1,num_frames);
%     for frIdx=1:num_frames;
%         try
%             L(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.bactLength*30;
%         catch
%             L(frIdx)=NaN;
%         end
%     end
%     for i=size(L,2):-1:1
%         if isnan(L(i))
%             L(i)=[];
%             T(i)=[];
%         end
%     end
%     [kL, l0, ~]  = fitLengthVsTime(T, L);
%     ElMat(cellIdx,:)=[kL, l0];
% end

% tic
% for analysisCycleId=1:9 % play with this
%     switch analysisCycleId
%         case 1
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160713 DC WT t0 30min\Results\together_05-12';
%         case 2
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160622 DC WT t0 22min\Results\together 3-5_9_10_12-15';
%         case 3
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160621 DC FtsZ ML159 35min\together 2_4-11_13_15';
%         case 4
%              PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160613 DC ML2159 t0 80min\Results\together 1_2_4-9';
%         case 5
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\161005 FtsW Mut 34min\Results\together 1-12\manual Tw 2';
%         case 6
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160816 WT DC 1305 FtsWGFP t0 50min\Results\together 01-06 and 10\manual Tw 2';
%         case 7
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160926 WT FtsW 18 488 37 min\Results\Together 01-07\manual Tw 2';
%         case 8
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160927 DC ML2159 FtsW t0 36min\Results\together_01-10\manual Tw 3';
%         case 9
%             PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160526 dc fOSFO T080\Results\together 1-9, 11,12,14-17';
%     end
%     load([PathNameResults '\tlMatrix.mat'])    
%     ElMat=getElMat(tlMatrix);    
%     save([PathNameResults '\ElMat.mat'], 'ElMat');
%     clear
% end
% duration=toc;
% fprintf('Finished! Elapsed time is:\n')
% breaktime(duration,'disp');