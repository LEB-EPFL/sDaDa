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

function [cellList, divededVar]= segCheckCellByCell(varargin)
% Input:
% single colour:
% (allSegmentedCellInfo, maskedImageCell, originalImNoMask, frIdx, divededVarPriv, PathNameResults, getMeshFig)
% dual colour:
% (allSegmentedCellInfo, maskedImageCell, originalImNoMask, frIdx, divededVarPriv, PathNameResults, getMeshFig, originalImageZ(:,:,frIdx),FtsZFig)
% 
narginchk(8,10);
allSegmentedCellInfo=varargin{1,1};
originalIm=varargin{1,2};
originalImNoMask=varargin{1,3};
frIdx=varargin{1,4};
divededVarPriv=varargin{1,5};
PathNameResults=varargin{1,6};
getMeshFig=varargin{1,7};
cellListPrev=varargin{1,8};
if size(varargin,2)==10;
    originalImageZ=varargin{1,9};
    FtsZFig=varargin{1,10};
    DualC=true;
else
    DualC=false;
end

    %% ALL CELL DETECTION/ MESH GENERATION PARAMETERS
    %image parameters
    param.areaMin = 80; 
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
    param.meshStep = 3; % 10
    param.meshTolerance = 0.01;
    param.meshWidth=45; % 50
    param.fsmooth=Inf;
    param.fmeshstep = 3; % 10

    % Limits on poition / number;
    param.roiBorder =5;
    param.noCellBorder = 5;
    param.imRoi = [1 1 Inf Inf];
    %%

    cellNum = size(allSegmentedCellInfo, 2)
    scrsz = get(0,'MonitorPositions');
    h1 = figure('Position',[scrsz(1,3)/2 50 scrsz(1,3)/2 scrsz(1,4)/2]);
    
    for cellIdx = 1: cellNum
        
        % If the bacteria exist in the raw image
        if ~isempty(allSegmentedCellInfo{1, cellIdx})
            % Create 'divededVarPriv' if it doesn't exist
            if isempty(divededVarPriv) || frIdx == 1
                divededVarPriv = zeros(cellNum, 3);
            end
            % If the cell is  NOT already divided analyse it
            if divededVarPriv(cellIdx, 3) == 0
                % The cell is not divided
                divided = 0;
                % Grab the area around the cell identified at the first step
                roiBoxTemp = allSegmentedCellInfo{1, cellIdx}.box;
                % Grab the contour of that cell
                contPriv = allSegmentedCellInfo{1, cellIdx}.contour;  % lastCont = segmentedCellInfo
                % Grab the centroid position
                centroidPrivTemp = allSegmentedCellInfo{1, cellIdx}.centroid;

                %% Expand the search ROI box by 40 pixels if possible
                % Check first if the cell is closer than 40 pixels to the
                % edge of the image. If so, make the roi smaller.
                if roiBoxTemp(1)-40 > 0 && roiBoxTemp(2)-40 > 0 && roiBoxTemp(1)+roiBoxTemp(3)+40 <= size(originalIm,2)&& roiBoxTemp(2) + roiBoxTemp(4)+40 <= size(originalIm,1)
                    roiBoxExp=[-40 -40 +80 +80];% [xmin ymin width height]
                    CentOffset=[40 40];
                else
                    roiBoxExp=[-40 -40 +80 +80];% [xmin ymin width height]
                    CentOffset=[40 40]; % [x y]
                    % if roiBoxTemp is closer than 40 px to border, match
                    % expansion to this distance.
                    if roiBoxTemp(1) < 40
                        roiBoxExp(1)=-roiBoxTemp(1);
                        CentOffset(1)=roiBoxTemp(1);
                    end
                    if roiBoxTemp(2) < 40
                        roiBoxExp(2)=-roiBoxTemp(2);
                        CentOffset(2)=roiBoxTemp(2);
                    end
                    if roiBoxTemp(1)+roiBoxTemp(3)+40 > size(originalIm,2)
                        roiBoxExp(3)=-roiBoxExp(1)+ size(originalIm,2)-(roiBoxTemp(1)+roiBoxTemp(3));
                    end
                    if roiBoxTemp(2)+roiBoxTemp(4)+40 > size(originalIm,1)
                        roiBoxExp(4)= - roiBoxExp(2)+ size(originalIm,1)-(roiBoxTemp(2)+roiBoxTemp(4));
                    end

                end
                
                centroidPriv = centroidPrivTemp + CentOffset;
                roiBox = roiBoxTemp + roiBoxExp;
                
                %% Crop ROI
                % Size and position of the crop rectangle, specified as a four-element position ...
                % vector of the form [xmin ymin width height].
                originalImCrop = imcrop(originalIm,roiBox);
                
                %% Resegment if the image is empty
                if sum(sum(originalImCrop))==0 
                    % redo the segmentation
                    originalImCrop = imcrop(originalImNoMask(:,:,frIdx),roiBox);
                    threshTmp=graythresh(originalImCrop);
                    maskTmp=im2bw(originalImCrop,threshTmp);
                    originalImCrop=originalImCrop.*maskTmp;
                end
                
                %% Redefine the new contour wrt the new ROI
                cont = rot90(contPriv,-2);
                cont(:, 1) = cont(:, 1) - roiBox(2) + 1;
                cont(:, 2) = cont(:, 2) - roiBox(1) + 1;
                
                %% Remove small regions
                % Identify individual blobs by seeing which pixels are connected to each other.
                % Each group of connected pixels will be given a label, a number, to identify it and distinguish it from the other blobs.
                % Do connected components labeling with either bwlabel() or bwconncomp().
                labeledImageCrop = bwlabel(originalImCrop, 8);     % Label each blob so we can make measurements of it

                % Get all the cells properties.  Can only pass in originalImage in version R2008a and later.
                cellMeasurements = regionprops(labeledImageCrop, originalImCrop, 'all');

                % Select certain cells based using the ismember() function.
                % Let's say that we wanted to find only those cells an area bigger 
                % than mean - sigma
                allCellAreas = [cellMeasurements.Area];

                % Get a list of the cells that meet our criteria and we need to keep.
                % These will be logical indices - lists of true or false depending on 
                % whether the feature meets the criteria or not.
                allowableAreaIndexes = allCellAreas >= 220; % Take the big objects.

                % Now let's get actual indexes, rather than logical indexes, of the features that meet the criteria.
                % keeperIndexes = find(allowableIntensityIndexes & allowableAreaIndexes);
                keeperIndexes = find(allowableAreaIndexes);

                % Extract only those blobs that meet our criteria, and
                % eliminate those cells that don't meet our criteria.
                % Note how we use ismember() to do this.  
                % Result will be an image - the same as labeledImage but with only the cells
                % listed in keeperIndexes in it.
                keeperCellsImage = ismember(labeledImageCrop, keeperIndexes);
                originalImCrop(~keeperCellsImage) = 0;  % Set all non-keeper pixels to zero.
                % find the centroid of the cells
                newCentroids = [cellMeasurements.Centroid];
                
                %% Select 1 region closest to previous position
                % If there is more than one cell in the FOV choose the one with the
                % center more close to the one in the previous frame
                if length(allCellAreas) > 1
                    % x y column
                    centroids = reshape( newCentroids, 2, size(newCentroids,2)/2)';

                    % (Xa - Xb)^2 and (Ya - Yb)^2  
                    distTemp = ( centroids - repmat(centroidPriv, size(centroids,1),1) ).*( centroids - repmat(centroidPriv, size(centroids,1),1) );
                    dist = sqrt(distTemp(:, 1) + distTemp(:, 2));
                    [~, minDistIdx] = min(dist);
                    centroid=centroids(minDistIdx,:);
                    if max(max(labeledImageCrop)) > 1
                        keeperCellsImage2 = ismember(labeledImageCrop, minDistIdx);
                        originalImCrop(~keeperCellsImage2) = 0;  % Set all non-keeper pixels to zero.
                    end
                else
                    centroid=reshape( newCentroids, 2, size(newCentroids,2)/2)';
                end
                
                %% Set-up for image classification
                % Dilation of maskedImageCell to be sure to include full region of the
                % cell (important for accurate FWHM).
                originalImCrop=expandCell(originalImCrop,originalImNoMask(:,:,frIdx),3, allSegmentedCellInfo{1,cellIdx}.box, roiBoxExp);
                
                res0 = originalImCrop; %original image
                resCur = originalImCrop;%current image
                resLast = originalImCrop;% undo/ backup image
                
                finished = false;
                succesfullyFinished=false;
                recordSteps=[];
                scrsz = get(0,'MonitorPositions');
                
                %% Try to skip user intervention if cell looks similar to previous frame
                % Check detected cell and skip the manual checks if all of the following are true:
                % - frIdx>1;
                % - There was an image of it in the previous time point.
                % exist('segmentedAllCellInfo{frIdx-1,1}.origCropImg','var')
                % - It was not too close to dividing in the previous time
                %   point. (Dmin/Dmax>threshold... 0.8)
                % - The detected image is not too different from the
                % previous one. (Hardest part)
                
                try
                    if isfield(cellListPrev{1,cellIdx},'origCropImg') % checks frIdx>1 and existence of image
                        if cellListPrev{1,cellIdx}.minD/cellListPrev{1,cellIdx}.maxD>0.8 % if not already constricting
                            % Now calculate difference
                            posDiff=centroid-cellListPrev{1,cellIdx}.centroid;
                            % [x y width height]
                            smallBox=[posDiff(1),... % x
                                posDiff(2),... % y
                                size(cellListPrev{1,cellIdx}.origCropImg,2)-1,... % width
                                size(cellListPrev{1,cellIdx}.origCropImg,1)-1]; % height
                            originalImCropSmall=imcrop(originalImCrop,smallBox);
                            
                            binPrev=cellListPrev{1,cellIdx}.origCropImg>0;
                            binNow=originalImCropSmall>0;
                            binDiff=xor(binPrev,binNow);
                            strEl=strel('disk',1);
                            binDiffEr=imerode(binDiff,strEl);
                            diffSum=sum(sum(binDiffEr));
                            areaPrev=sum(sum(binPrev));
                            areaNow=sum(sum(binNow));
                            
                            % if (eroded) per pxl difference is less than 1.5% (check
                            % these values)
                            if diffSum/areaPrev<0.015
                                finished=true;
                            else
                                finished=false;
                            end
                        end
                    end
                catch errorInAutoSeg
                    finished=false;
                    disp(['segCheckCellByCell line 283 ' errorInAutoSeg.identifier]);
                end
                
                %%
                
                if ~finished % open fig if segmentation not finished
                    h2 = figure('Position',[8 scrsz(1,4)*0.041 scrsz(1,3)/2 scrsz(1,3)*0.45]); % Figure for plotting cell image and contour and deciding what to do with it (split,delete,divided,...)
                    title('h2');
                    fprintf('\n\n');
                end
                
                while ~succesfullyFinished %
                while ~finished % User-aided segmentation

                  %% update fig
                  figure(h2);
                  imagesc(resCur);
                  colormap(jet)
                  axis equal,
                  hold on, 
                  plot(cont(:,2), cont(:,1), 'r--', 'Linewidth', 2)
                  
                  format long g;
                  format compact;
                  textFontSize = 10;	% Used to control size of "cell number" labels put atop the image.
                  labelShiftX = -7;	% Used to align the labels in the centers of the cell.

                  % Put the labels on the rgb labeled image
                  if length(allCellAreas) > 1
                    text(centroids(minDistIdx,1) + labelShiftX, centroids(minDistIdx,2), num2str(cellIdx), 'FontSize', textFontSize, 'FontWeight', 'Light','Color','w' );
                  end
                  text(centroidPriv(1,1) + labelShiftX, centroidPriv(1,2), num2str(cellIdx), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','w' );
                  %% Prompt user to pick an option
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
                  fprintf(' (10) Cell Divided\n');
                  fprintf(' (0) Expand cell\n');
                  figure(h2); % make figure h2 active
                  str = input('Choose 0-10: ','s');
                  recordSteps = [recordSteps,str2num(str)]
                  switch str
                  case '1' %Split cells
                    resLast = resCur;
                    try
                    resCur = splitCell(resCur,param.areaMin,h2);
                    catch
                        disp('Whoopsie, something went wrong there! Try again, but harder!')
                    end
                  case '2' % Join Cells
                    resLast = resCur;
                    try
                    resCur = joinCell(resCur,h2);
                    catch
                        disp('Whoopsie, I am not sure what happened there! Could you try that again?')
                    end
                  case '3'
                    resLast = resCur;
                    try
                    resCur = deleteCell(resCur,h2);
                    % divided = 2;
                    catch
                        disp('Not today...')
                    end
                  case '4'
                    fprintf('\n1');
                    resLast = resCur;
                    try
                    resCur = selectGoodCell(resCur,h2);
                    catch
                        disp('This is not the cell you are looking for!')
                    end
                  case '5'
                    resLast = resCur;
                    try
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
                    catch
                        disp('I am afraid I cannot do that, Dave.')
                    end
                  case '6'
                    resLast = resCur;
                    try
                    resCur = imclearborder(resCur);
                    catch
                        disp('You have crossed the line, pall.')
                    end
                  case '7'
                    fprintf('Undid last change\n');
                    resCur = resLast;
                    try
                        if recordSteps(end-1) == 10
                            divided=0;
                        end
                    catch
                    end
                  case '8'
                    fprintf('Undid all changes\n');
                    resCur = res0;
                    divided=0;
                  case '9'
                    %check if there is only one "object" selected
                    TestCur=bwlabel(resCur);
                    if max(max(TestCur)) < 2
                        finished = true;
                    else
                            disp('There are multiple objects in the FOV, pick one')
                    end
                  case '10'
                     try
                         resLast = resCur;
                         resCur = deleteCell(resCur,h2);
                         divided = 1;
                     catch
                         disp('Failed dividing the cells...')
                     end
                  case '0'
                         try
                             resLast = resCur;
                             fprintf('\nEnter radius of strel for smoothing via image closing\n');
                             fprintf('Default is 3, press Enter for default\n');
                             str = input('Strel radius: ','s');
                             if strcmp(str,'')
                                closeNum = 3;
                             else
                                closeNum = str2num(str);
                             end
                             resCur=expandCell(resCur,originalImNoMask(:,:,frIdx),closeNum, allSegmentedCellInfo{1,cellIdx}.box, roiBoxExp);
                         catch
                             disp('Expansion failed.')
                         end
                  otherwise
                    fprintf('\nInvalid choice\n\n');
                  end
                end % end sub-segmentation
                %% Segmenattion supposedly finished, proceed to measurement
                cellIm = resCur;
                if sum(sum(cellIm),2)==0 % Nothing in resCur/cellIm/h2
                    %% Cell is deleted or divided
                    try
                        close(h2)
                    catch
                    end
                     disp('frame:')
                     disp(frIdx)
                     cellList{1, cellIdx} = [];
                     succesfullyFinished=true;
                else
                    %% Measurements of MTS (and divisome) channel
                     disp('frame:')
                     disp(frIdx)
                     try
                         try % Update originalImcrop based on resCur, for getMesh...
                             % resCur can be logical or  masked intensity
                             % image. (e.g. after smoothing)
                             originalImCrop=expandCell(resCur,originalImNoMask(:,:,frIdx),0, allSegmentedCellInfo{1,cellIdx}.box, roiBoxExp);
                         catch
                             warning('Failed to update originalImCrop based on resCur.')
                         end
                         % Measurement on MTS channel
                         cellList{1, cellIdx} = getMeshAlgForSIM(originalImCrop, cellIm, param, cellIdx, frIdx, PathNameResults, getMeshFig);
                         if DualC % If DualC, plot Z channel, ... and measure Z-ring
                             %% Based on MTS, define ROI and Mask to inspect in divisome channel
                             % create mask based on cell image, zero becomes false, other values true.
                             cellMask=cellList{1,cellIdx}.origCropImg>0;
                             % expand mask, to include 
                             se=strel('disk',5);
                             cellMask=imdilate(cellMask,se);
                             % expand "mask borders" such that if fits the size
                             % roiBox. (+1)
                             tmp=zeros(floor(cellList{1,cellIdx}.box(4)+roiBoxExp(4)+1), floor(cellList{1,cellIdx}.box(3)+roiBoxExp(3)+1)); % this doesnt work when cell is close to border: roiBox(4)+1,roiBox(3)+1
                             y1=-roiBoxExp(1)+1;
                             x1=-roiBoxExp(2)+1;
                             y2=size(cellMask,2)+y1-1;
                             x2=size(cellMask,1)+x1-1;
                             tmp(x1:x2,y1:y2)=cellMask;
                             cellMask=tmp;
                             try % adjust roiBoxZ and CellMask and mask FtsZ image
                                 roiBoxZ=cellList{1,cellIdx}.box;
                                 % now expand roiBoxZ
                                 roiBoxZ=roiBoxZ+roiBoxExp;
                                 % set the x and y coordinates of roiBoxZ
                                 % relative to image instead of relative to
                                 % roiBox.
                                 roiBoxZ(1)=roiBoxZ(1)+roiBox(1);
                                 roiBoxZ(2)=roiBoxZ(2)+roiBox(2);
                                 originalImZCrop=imcrop(originalImageZ,roiBoxZ);
                                 cellMask=imcrop(cellMask,[0 0 size(originalImZCrop,2) size(originalImZCrop,1)]);
                                 originalImZMasked=originalImZCrop.*cellMask;
                             catch
                                 disp('Not this error again... originalImZMasked=originalImZCrop.*cellMask')
                                 %% This part fixes an error when the cell
                                 % touches the edge of the roiBox.
                                 try
                                     originalImZMasked=originalImZCrop(1:end-1,1:end-1).*cellMask;
                                 catch
                                     try
                                         originalImZMasked=originalImZCrop(1:end-1,:).*cellMask;
                                     catch
                                         try
                                             originalImZMasked=originalImZCrop(:,1:end-1).*cellMask;
                                         catch
                                             disp('Even this did not work.')
                                         end
                                     end
                                 end
                                 
                             end
                             roiCor=[-cellList{1,cellIdx}.box(1)-roiBoxExp(1),-cellList{1,cellIdx}.box(2)-roiBoxExp(2)];
                             
                             %% Select/create figure to plot divisome measurements
                             try
                                 figure(FtsZFig);
                             catch
                                 FtsZFig=figure;
                             end
                             figure(FtsZFig),
                             
                             %% get cLine from MTS for divisome channel measurements
                             cLine=getCLine(cellList{1,cellIdx});
                             cLine(:,1)=cLine(:,1)+roiCor(1); % correct position relatively to the roi box
                             cLine(:,2)=cLine(:,2)+roiCor(2);
                             
                             %% Find and measure the divisome
                             try
                                cellList{1,cellIdx} = measureZ(originalImZMasked, cellList{1,cellIdx}, cLine, roiCor, cellIdx, frIdx, PathNameResults,FtsZFig);
                             catch ME
                                 disp('u got divispwned')
                                 fprintf([ME.identifier '\n'])
                             end
                         end
                         
                         try
                            close(h2)
                         catch
                         end
                         
                         succesfullyFinished=true;
                     catch ME % If the measurements, go back to segmentation
                         disp('OMGZ something went wrong!')
                         fprintf([ME.identifier '\n'])
                         finished=false;
                         succesfullyFinished=false;
                         try % reselect/create figure
                            figure(h2);
                         catch
                             h2=figure('Position',[8 scrsz(1,4)*0.041 scrsz(1,3)/2 scrsz(1,3)*0.45]);
                         end
                     end
                end
                end %end "succesfullyFinished"
                %% Update divededVar
                divededVar(cellIdx, :) = [frIdx cellIdx divided]; 
                
                %% Debug: Save variables for testing autosemiauto
                try
                    cellList{1,cellIdx}.diffSum=diffSum;
                    cellList{1,cellIdx}.areaPrev=areaPrev;
                    cellList{1,cellIdx}.areaNow=areaNow;
                    cellList{1,cellIdx}.recordSteps=recordSteps;
                catch
                end
                
            else % if is already divided remember it
                divededVar(cellIdx, :) = [frIdx cellIdx 1]; 
                cellList{1, cellIdx} = []; 
            end  % end if the cell is already divided
                    
        else % end if the cell does not exist
            cellList{1, cellIdx} = [];
        end
          
    end % end cell cycle
    try
        close(h1)
    catch
        disp 'cannot close figure h1, does it even exist?'; %I added this because closing h1 causes errors when it doesn't exist
    end
end

%----------------------------------------
function res = splitCell(res,minCellSize,h2)
figure(h2);
res = watershedSeedSeg4(res,minCellSize);
end
%----------------------------------------
function res0= joinCell(res,h2)
figure(h2);
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

joinCellIm= logical(zeros(size(res)));
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
figure(h2);
[x,y] = ginputc('ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res,4); % Connectivity 4 to remove diagonally connected pixels
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
zeroIm = false(size(res));
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

function cLine=getCLine(cellList)
mesh=cellList.mesh;
cLinetmp = zeros(size(mesh,1),2);
cLinetmp(:,1) = mean([mesh(:,1),mesh(:,3)],2);
cLinetmp(:,2) = mean([mesh(:,2),mesh(:,4)],2);
% smooth the centreLine
[x,y]=spline2d(cLinetmp(:,1),cLinetmp(:,2),3); %nbreak=3
cLine=[x(:),y(:)];
end

function [xOut,yOut]=spline2d(x,y,nBreak)
nPt = numel(x);
nPtOut = 10*nPt;% upsample the output, but keep discrete

t=linspace(0,1,nPt);
tOut=linspace(0,1,nPtOut);
ppX = splinefit(t,x,nBreak);
ppY = splinefit(t,y,nBreak);
xOut = ppval(ppX,tOut);
yOut = ppval(ppY,tOut);
end

function expCur=expandCell(resCur,originalImNoMask,radius, roiBoxSmall, roiBoxExp)
% Function to expand the roi/contour of a cell. 

% Turn image into mask
mask=resCur>0;

% Expand mask
se=strel('disk',radius);
mask=imdilate(mask,se);

% Adjust roiBoxBig and mask, then mask image
roiBoxBig=roiBoxSmall+roiBoxExp;
% now expand roiBoxBig
% roiBoxBig=roiBoxBig+roiBoxExp;
% set the x and y coordinates of roiBox
% relative to image instead of relative to
% roiBox.
roiBoxBig(1)=roiBoxBig(1); %+roiBoxSmall(1);
roiBoxBig(2)=roiBoxBig(2); %+roiBoxSmall(2);
originalImCrop=imcrop(originalImNoMask,roiBoxBig);
mask=imcrop(mask,[0 0 size(originalImCrop,2) size(originalImCrop,1)]); %just to be sure the sizes match
expCur=originalImCrop.*mask;
end
