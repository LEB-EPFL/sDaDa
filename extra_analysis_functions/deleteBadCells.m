function [cellInfo, divededVarEachFrm] = deleteBadCells(cellInfo, divededVarEachFrm, badCells)
% This is a function to delete bad cells. Use with caution.
% Workflow: first manually look for bad cells, figure out their number
% (cellInfo{frame,1}{1,cell number})
% Input their numbers in the badCells vector.
% run section "Checking" to plot WvT and LvT for these cells to make sure
% you want to delete them.
% The actual "Deleting", deletes cells from cellInfo and DivededVarEachFrm,
% which results in new numbers for the cells. It is as if they were never
% there!
% "Double checking": plot the remaining cells.

%% Specify bad cells
%badCells = [8 30]; % Input
badCells=reshape(badCells,[1,numel(badCells)]);
badCells = flip(sort(badCells)); % sort the numbers in descending order.
% It is important to delete cells in descending order, since deleting a
% cell will decrease the number of all following cells.
%% Checking
% Took code from "PlotFit10Cells"
T=0:5:(75-1)*5;
figure,
hold on
for cellIdx=badCells;
    W=zeros(1,75); % preallocation
    for frIdx=1:75;
        try
            W(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
        catch
            W(frIdx)=0;
        end
    end
    plot(T,W,'.:','MarkerSize',15);
    [~, ~, ~, fitOutTimeNoBinNoLengF]  = fitFeingoldModBound(T, W);
    plot(fitOutTimeNoBinNoLengF(:, 1), fitOutTimeNoBinNoLengF(:, 2), '-') % fit with Feingold
end
hold off

figure,
hold all
for cellIdx=badCells;
    W=zeros(1,75); % preallocation
    L=zeros(1,75);
    for frIdx=1:75;
        try
            W(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
            L(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.bactLength*30;
        catch
            W(frIdx)=0;
            try
                if W(frIdx-1)==0;
                    L(frIdx)=NaN;
                else
                    L(frIdx)=L(frIdx-1);
                end
            catch
                L(frIdx)=NaN;
            end
        end
    end
    plot(L,W,'.:','MarkerSize',15); % ,'.','MarkerSize',15
    
end
title('Normalized waist width vs Length. Bad cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

ctu=input('Continue? y/n \n');
if ~strcmp(ctu,'y')
    return
end

%% Deleting
numFrames=size(cellInfo,1);
for ii = 1:length(badCells)
    cellIdx=badCells(ii);
    % delete cell from cellInfo
    for frIdx=1:numFrames
        cellInfo{frIdx,1}(1,cellIdx:end-1)=cellInfo{frIdx,1}(1,cellIdx+1:end); % move all fields forward, overwriting the bad cell
        cellInfo{frIdx,1}=cellInfo{frIdx,1}(1,1:end-1); % clear last field (it is a duplicate now)
    end
    
    % delete cell from divededVarEachFrm
    cellDVIdx=find(divededVarEachFrm(:,2) == cellIdx); % find the cell
    cellDVIdx=flip(sort(cellDVIdx));
    for jj= cellDVIdx
        divededVarEachFrm(jj,:)=[];
        laterCellIdx=find(divededVarEachFrm(:,2) > cellIdx); % find the cells that need their cell number decreased
        divededVarEachFrm(laterCellIdx,2)=divededVarEachFrm(laterCellIdx,2)-1;
    end
end
disp('Cells deleted!')
%% Double Checking
% Took code from "PlotFit10Cells"
T=0:5:(75-1)*5;
figure,
hold on
for cellIdx=1:size(cellInfo{1,1},2);
    W=zeros(1,75); % preallocation
    noplot=false;
    for frIdx=1:75;
        try
            W(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
            WMat(cellIdx,frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
        catch
            try
                W(frIdx)=0;
                WMat(cellIdx,frIdx)=0;
            catch
                fprintf('Cell %i does not exist', cellIdx)
                noplot=true;
            end
        end
    end
    if ~noplot
        plot(T,W,'.:','MarkerSize',15);
    end
end
hold off

figure,
hold all
for cellIdx=1:size(cellInfo{1,1},2);
    W=zeros(1,75); % preallocation
    L=zeros(1,75);
    for frIdx=1:75;
        try
            W(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
            L(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.bactLength*30;
        catch
            W(frIdx)=0;
            if W(frIdx-1)==0;
                L(frIdx)=NaN;
            else
                L(frIdx)=L(frIdx-1);
            end
        end
    end
    plot(L,W,'.:','MarkerSize',15);
end
title('Normalized waist width vs Length. 10 cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

ctu=input('Continue? y/n \n');
if ~strcmp(ctu,'y')
    return
end
