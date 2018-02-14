% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Anna Archetti and Aster Vanhecke 
% 
% Contact:
% e-mail:
% anna.archetti@epfl.ch
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

function cellList = semiAutoSegForSIM(originalIm, PathNameResults, getMeshFig)
% function semiAutoSegForSIM, to perform segmentation on the first frame.
% format: cellList = semiAutoSegForSIM(originalIm, PathNameResults, getMeshFig)

    %% ALL CELL DETECTION/ MESH GENERATION PARAMETERS
    %image parameters
    param.areaMin = 1000; 
    param.closeNum = 10;
    
    param.invertimage=1;
    param.erodeNum = 0;
    %segmentation parameters
    param.edgemode='sauvola';
    param.localNeighbourhood=[15 15];
    param.localThresh = 0.2;
    param.smoothNum= 10;
    

    % Shape parameters
    param.interpoutline =0;
    param.interpWeights =[0.5 0.5];
    param.interpSigma =[5 5];
    param.thresFactorF = 1;

    % Mesh parameters;
    param.getmesh = 1;
    param.meshStep = 3;% 10
    param.meshTolerance = 0.01;
    param.meshWidth=45;% 50
    param.fsmooth=Inf;
    param.fmeshstep = 3;%10

    % Limits on poition / number;
    param.roiBorder =5;
    param.noCellBorder = 5;
    param.imRoi = [1 1 Inf Inf];
    %%
    
    res0= originalIm; %original image
    resCur = originalIm;%current image
    resLast = originalIm;% undo/ backup image

    finished = false;
    recordSteps=[];
    scrsz = get(0,'MonitorPositions');
    h2 = figure('Position',[0 scrsz(1,4)/10 scrsz(1,3)/2 scrsz(1,3)/2]);
    title('h2')
    fprintf('\n\n');
    while ~finished
      %update fig
      figure(h2);
      labelIm = bwlabel(resCur);
      RGB  = label2rgb(labelIm,'jet','k','shuffle');
      imshow(RGB,'InitialMagnification',60);

      %options
      fprintf('Cell segmentation options:\n');
      fprintf(' (1) Split cells\n');
      fprintf(' (2) Join cells\n');
      fprintf(' (3) Delete cells\n');
      fprintf(' (4) Select good cells (deletes other cells)\n');
      fprintf(' (5) Smooth cells\n');
      fprintf(' (6) Delete edge cells\n');
      fprintf(' (7) Undo last change \n');
      fprintf(' (8) Discard all changes (reset to original image\n');
      fprintf(' (9) Finished segmentation\n');
      str = input('Choose 1-9: ','s');
      recordSteps = [recordSteps,str2num(str)]
      switch str
      case '1' %Split cells
        resLast = resCur;
        resCur = splitCell(resCur,param.areaMin,h2);
      case '2' % Join Cells
        resLast = resCur;
        resCur = joinCell(resCur,h2);
      case '3' % Delete Cell
        resLast = resCur;
        resCur = deleteCell(resCur,h2);
      case '4' % Select good cells
        fprintf('\n1');
        resLast = resCur;
        resCur = selectGoodCell(resCur,h2);
      case '5' % Smooth cells
        resLast = resCur;
        fprintf('\nEnter radius of strel for smoothing via image closing\n');
        fprintf('Default is %d, press Enter for default\n',param.closeNum);
        str = input('Strel radius: ','s');
        if strcmp(str,'')
          closeNum = param.closeNum;
        else
          closeNum = str2num(str);
        end
        resCur= smoothCell(resCur,closeNum);
        recordSteps = [recordSteps,closeNum];
      case '6' % Clear cells touching border
        resLast = resCur;
        resCur = imclearborder(resCur);
      case '7' % Undo last
        fprintf('Undid last change\n');
        resCur = resLast;
      case '8' % Reset/Undo all
        fprintf('Undid all changes\n');
        resCur = res0;
      case '9' % Finish
        cellIm = resCur;
        finished = true;
      otherwise
        fprintf('\nInvalid choice\n\n');
      end
    end
    mkdir(PathNameResults);
    print(h2,'-dpng',[PathNameResults '\AllCellsSegmented.png']);
    cellList = getMeshAlgForSIM(originalIm, cellIm, param, 0, 0, PathNameResults, getMeshFig);
end

%----------------------------------------
function res = splitCell(res,minCellSize,h2)

 figure(h2);
res = watershedSeedSeg4(res,minCellSize);
end
%----------------------------------------
function res0= joinCell(res,h2)
figure(h2)
[x,y] = ginputc(2,'ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

%get the selected cells
labelIm = bwlabel(res);
nCell = numel(xI);
joinCell = zeros(nCell,1);
for ii = 1:nCell;
  labelIm = bwlabel(res);
  joinCell(ii) = labelIm(yI(ii),xI(ii));
end

joinCellIm= false(size(res));
joinCell = unique(joinCell);
for ii = 1:numel(joinCell)
  joinCellIm(find(labelIm==joinCell(ii)))=1;
end

%close the image until they are connected
isJoined=false;
ii=1;
while ~isJoined
  joinCellLabel = bwlabel(joinCellIm);
  if max(joinCellLabel(:))<=1
    isJoined=true;
  else
    joinCellIm = imclose(joinCellIm,strel('disk',ii));
    ii=ii+1;
  end
end

%copy the joined image back to the main image
res0=res;
res0(joinCellIm==1)=1;
end

%----------------------------------------
function goodCellIm= selectGoodCell(res,h2)
figure(h2)
[x,y] = ginputc('ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res);
nCell = numel(xI);
goodCell = zeros(nCell,1);
for ii = 1:nCell;
  labelIm = bwlabel(res);
  goodCell(ii) = labelIm(yI(ii),xI(ii));
end

goodCellIm= false(size(res));
goodCell = unique(goodCell);
for ii = 1:numel(goodCell)
  goodCellIm(find(labelIm==goodCell(ii)))=1;
end
end

%----------------------------------------
function smoothCellIm= smoothCell(res,closeNum)

labelIm = bwlabel(res);
zeroIm = logical(zeros(size(res)));
oneCellIm = zeroIm;
smoothCellIm=zeroIm;

nCell = max(labelIm(:));
for ii = 1:nCell
  oneCellIm = zeroIm;
  oneCellIm(labelIm==ii)=1;
  %smooth
  oneCellIm= imclose(oneCellIm,strel('disk',closeNum));
  smoothCellIm(oneCellIm==1)=1;
end
end

%----------------------------------------
function badCellIm= deleteCell(res,h2)
figure(h2)
[x,y] = ginputc('ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res);
nCell = numel(xI);
badCell = zeros(nCell,1);
for ii = 1:nCell
  labelIm = bwlabel(res);
  badCell(ii) = labelIm(yI(ii),xI(ii));
end

badCellIm= res;
badCell = unique(badCell);
for ii = 1:numel(badCell)
  badCellIm(find(labelIm==badCell(ii)))=0;
end
end

