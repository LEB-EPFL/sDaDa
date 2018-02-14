function cellList = getMeshAlgForSIM(img, regions0, p, cellIdx, frIdx, PathNameResults,getMeshFig)
    % cellList = getMeshAlgForSIM(img, regions0, p, cellIdx, frIdx, PathNameResults,getMesh)
    % getMeshFig is the figure handle for the figure where the measured
    % diameter, the contour, mesh etc... are plotted.
    % Based on getMeshAlg from MicrobeTracker
    % 
    % Authors: Anna Archetti and Aster Vanhecke
    
    %% Set-up parameters etc
    format long g;
    format compact;
    
    warning('off', 'arguments:nany') % supresses warningspam from splinefit in getFWXM
    
    % Initialize strel
    [se, ~, ~] = initP();

    imsizes=[size(img,1) size(img,2)];

    regions0 = bwlabel(regions0,4);

    stat0 = regionprops(regions0,'area','boundingbox');
    cell0 = 0;
    frame = 1;
    regmax = length(stat0);
    
    loopVec=1:regmax;
    
    %% Check for diagonally connected pixels by removing small regions
    for xoxo=1:regmax
        if stat0(xoxo).Area<5 % small area
            loopVec(xoxo)=[];
        end
    end
    
    %% Function body
    for reg=loopVec % Is usually just one region.
        %% Initialize ROI, mask etc. for current reg
        statC = stat0(reg);

        % otherwise compute the properties for proper splitting
        roiBox(1:2) = max(statC.BoundingBox(1:2)-p.roiBorder,1); % coordinates of the ROI box
        roiBox(3:4) = min(statC.BoundingBox(1:2)+statC.BoundingBox(3:4)+p.roiBorder,[size(img,2) size(img,1)])-roiBox(1:2);
        roiRegs = imcrop(regions0,roiBox); % ROI with all regions labeled

        roiMask = roiRegs==reg; % ROI with current region #reg labeled
        roiLabeled = roiRegs.*roiMask;
        roiImg = imcrop(img,roiBox);
        roiOrigIm = roiImg.*roiMask;

        % Get all the cell properties.  
        cellMeasurements = regionprops(roiLabeled, roiOrigIm, 'all');
        
        numberOfCells = size(cellMeasurements, 1);

        %initialize cellList
        cellList{frame, numberOfCells} = cell(1);

        % Get and smooth the boundary of the region
        [ii,jj]=find(bwperim(roiMask),1,'first');
        cCell=bwtraceboundary(roiMask,[ii,jj],'n',8,inf,'counterclockwise');
        
        % Interpolate outline
        if isfield(p,'interpoutline') && p.interpoutline
          ctrlist = interpoutline(cCell,roiImg,p);
        else
          ctrlist = {cCell};
        end

        if isempty(ctrlist), display(['fitting region ' num2str(reg) ' - failed, no outline creating during smoothing']); end
       %% 
        for cind = 1:length(ctrlist)

            cCell = ctrlist{cind};
            if p.fsmooth>0 && p.fsmooth<Inf
                fpp = frdescp(cCell);
                cCell = ifdescp(fpp,p.fsmooth);
                cCell=[cCell;cCell(1,:)]; % make the last point the same as the first one
            end
            % Rotation and translation
            cCell = rot90(cCell,2);
            cCell(:,1) = cCell(:,1)+roiBox(1)-1;
            cCell(:,2) = cCell(:,2)+roiBox(2)-1;

            %% MESH DEFINITION  
            % Save the data
            if p.getmesh, mesh = model2mesh(cCell,p.meshStep,p.meshTolerance,p.meshWidth); end

            if (~p.getmesh || length(mesh)>4) && ...
             min(cCell(:,1))>= p.imRoi(1) &&...
             min(cCell(:,2))>= p.imRoi(2) &&...
             max(cCell(:,1))<=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)) && ...
             max(cCell(:,2))<=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)) 

                cellStruct.algorithm = 1;
                if p.getmesh, cellStruct.mesh = mesh; end

                %% Find the furthest away points(the cell poles)
                numberOfCoords = size(cCell,1);
                maxDistance = zeros(1,numberOfCoords); % preallocation
                indexOfMax = zeros(1,numberOfCoords);
                x = cCell(:,2);
                y = cCell(:,1);
                for k = 1 : numberOfCoords
                    distances = sqrt((x-x(k)).^2 + (y-y(k)).^2);
                    [maxDistance(k), indexOfMax(k)] = max(distances);
                end
                pole1PosIdx = find(maxDistance == max(maxDistance));
                pole2PosIdx = indexOfMax(pole1PosIdx);

                %% compute and smooth centreLine
                nBreak = 3;
                cLine = zeros(size(mesh,1),2);
                cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
                cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);
                % smooth the centreLine
                [x,y]=spline2d(cLine(:,1),cLine(:,2),nBreak);

                cLineSmooth=[x(:),y(:)];

                %% Rotate mesh and contour for plotting
                cont = rot90(cCell,-2);
                cont(:, 1) = cCell(:, 1)-roiBox(1)+1;
                cont(:, 2) = cCell(:, 2)-roiBox(2)+1;

                mesh1 = rot90(mesh(:,1:2), -2);
                mesh2 = rot90(mesh(:,3:4), -2);
                mesh1(:,1) = mesh1(:,1)-roiBox(2)+1;
                mesh1(:,2) = mesh1(:,2)-roiBox(1)+1;
                mesh2(:,1) = mesh2(:,1)-roiBox(2)+1;
                mesh2(:,2) = mesh2(:,2)-roiBox(1)+1;

                mesh1=flipdim(mesh1,2);
                mesh2=flipdim(mesh2,2);
                cellCentroid = cellMeasurements(reg, 1).Centroid;
                
                %% Plot contour, mesh, poles and central line on masked image.
                try 
                    getMeshFig.BeingDeleted; % Test whether the figure is deleted. (Cannot find a more sensible way to do it :s )
                catch
                    getMeshFig=figure;
                end
                figure(getMeshFig);
                hold off,
                subplot(1,2,1), 
                imshow(roiOrigIm), hold on,
                plot(cont(:,1),cont(:,2), 'b--'), 
                plot(mesh1(:,1),mesh1(:,2),'k'),
                plot(mesh2(:,1),mesh2(:,2),'k'),
                plot( cCell(pole1PosIdx,1)-roiBox(1)+1, cCell(pole1PosIdx,2)-roiBox(2)+1,'ro', ...
                      cCell(pole2PosIdx,1)-roiBox(1)+1, cCell(pole2PosIdx,2)-roiBox(2)+1,'ro'),
                plot( cLineSmooth(:,1)-roiBox(1)+1, cLineSmooth(:,2)-roiBox(2)+1, 'r--')
                format long g;
                format compact;
                
                textFontSize = 10;	% Used to control size of "cell number" labels put atop the image.
                labelShiftX = -7;	% Used to align the labels in the centers of the cell.
                % Put the labels on the rgb labeled image also.
                if cellIdx == 0
                    text(cellCentroid(1,1) + labelShiftX, cellCentroid(1,2), num2str(reg), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','w' );
                else
                    text(cellCentroid(1,1) + labelShiftX, cellCentroid(1,2), num2str(cellIdx), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','w' );
                end
                xlabel('[pixel]')
                ylabel('[pixel]')
                axis equal
                legend('contour', 'mesh','mesh', 'pole', 'central line') 
                
                %% Compute length of the bacterium
                [bactLength, seglen]= arclength(x, y, 'spline');
                len=[0; cumsum(seglen)];
                
                %% PARAMETER definition
                % find points at which to crop the data and plot the histograms
                lStep = 1;
                lWidth = 3; % 1 creates pixel size-related artefacts. 3 Works.
                l0 = 0:lStep:max(len);
                lMid = l0(1:end-1)+diff(l0)/2;
                nStep = numel(l0)-1;    
                lengthOut = lMid';   
                histStep = 0.5;

                % Variables initialization---------------------------------
                diamFWHM = [];
                diamMesh = [];
                meshLim = [];
                xy1ForIm = [];
                xy2ForIm = [];
                %preallocation
                intProfCelltmp{2,nStep}=[];
                c0AllPos(nStep, 2)=0;
                meshLim(nStep,2) =0;
                diamMesh(nStep) =0;
                
                %% Step through points on centreLine and measure width
                for ii = 1:nStep

                    % Find current mesh centre position, and +- lWidth/2
                    % Probably best to smooth this!
                    cPt(:,1) = interp1(len,cLineSmooth(:,1),[l0(ii), l0(ii)+lWidth/2, l0(ii)+lWidth]');
                    cPt(:,2) = interp1(len,cLineSmooth(:,2),[l0(ii), l0(ii)+lWidth/2, l0(ii)+lWidth]');
                    c0Temp = cPt(2,:);%this is the midpoint
                    
                    % Save the C0tot posisition of all the central position
                    % along the bacteria length during slicing
                    c0AllPos(ii, :) = c0Temp;

                    %% get angle for rotation
                    dX = cPt(3,1)-cPt(1,1);
                    dY = cPt(3,2)-cPt(1,2);

                    %  Get meshMin Max distances ------------------------
                    % meshLim = [-meshBelow, meshAbove];
                    meshLim(ii,:) = getCellDist(mesh,c0Temp,dY,dX);
                    diamMesh(ii) = meshLim(ii, 2) - meshLim(ii, 1);
                    % ----------------------------------------------------
                    
                    theta = atan2(dY,dX);

                    if ~any(isnan(meshLim(ii, :)))
                        
                        rowNum = size(roiOrigIm, 1);
                        colNum = size(roiOrigIm, 2);
                        
                        c0(1)= c0Temp(1)-roiBox(1)+1;
                        c0(2)= rowNum - (c0Temp(2)-roiBox(2)+1);
                        
                        distYPointFirstRow = rowNum - c0(2);                       
                        distXPointLastCol = colNum - c0(1);
                         
                        %% Image translation + rotation
                        % Such that current point on cLine becomes center,
                        % and image is aligned to cLine in this point.
                        % NOTE I changed floor in round my function trans
                        if distYPointFirstRow*2 >= rowNum && c0(1)*2 >= colNum
                            imCent = zeros(round(distYPointFirstRow*2)+1, round(c0(1)*2)-1);
                            imCent(1: rowNum, 1:colNum) = roiOrigIm;
                        end
                        
                        if distYPointFirstRow*2 < rowNum && c0(1)*2 < colNum
                            if round( c0(2) )*2 - 1 > rowNum 
                                imCent = zeros(round(c0(2)*2)-1, round(distXPointLastCol*2)+1);
                                newColNum = round(distXPointLastCol*2)+1;
                                newRowNum = round( c0(2) )*2 - 1;
                                imCent(newRowNum - rowNum + 1 : newRowNum - rowNum + 1 + rowNum -1, ...
                                newColNum - colNum + 1 : newColNum - colNum + 1 + colNum-1) = roiOrigIm;
                            else
                                imCent = zeros(round(c0(2)*2)+1, round(distXPointLastCol*2)+1);
                                newColNum = round(distXPointLastCol*2)+1;
                                newRowNum = round( c0(2) )*2 +1;
                                imCent(newRowNum - rowNum + 1 : newRowNum - rowNum + 1 + rowNum -1, ...
                                newColNum - colNum + 1 : newColNum - colNum + 1 + colNum-1) = roiOrigIm;
                            end
                        end
                        
                        if distYPointFirstRow*2 < rowNum && c0(1)*2 >= colNum
                            if round( c0(2) )*2 - 1 > rowNum 
                                imCent = zeros(round(c0(2)*2)-1, round(c0(1)*2)-1);
                                newRowNum = round(c0(2)*2)-1;
                                imCent(newRowNum - rowNum + 1 : newRowNum - rowNum + 1 + rowNum -1, 1:colNum) = roiOrigIm;
                            else
                                imCent = zeros(round(c0(2)*2)+1, round(c0(1)*2)-1);
                                newRowNum = round(c0(2)*2)+1;
                                imCent(newRowNum - rowNum + 1 : newRowNum - rowNum + 1 + rowNum -1, 1:colNum) = roiOrigIm;
                            end
                        end
                        
                        if distYPointFirstRow*2 >= rowNum && c0(1)*2 < colNum
                            imCent = zeros(round(distYPointFirstRow*2)+1, round(distXPointLastCol*2)+1);
                            newColNum = round(distXPointLastCol*2)+1;
                            
                            imCent(1: rowNum, newColNum - colNum + 1 : newColNum - colNum + 1 + colNum-1) = roiOrigIm;
                        end
                        
                        % Rotate parallel to xAxis 
                        thetaRot = -theta;
                        R = [cos(thetaRot),-sin(thetaRot );...
                        sin(thetaRot ),cos(thetaRot)];
                    
                    % Mesh Translation to center in c0
                    xy1Temp = [mesh(:, 1) - roiBox(1) + 1 - c0(1),   mesh(:, 2) - roiBox(2) + 1 + c0(2) - rowNum];
                    xy2Temp = [mesh(:, 3) - roiBox(1) + 1 - c0(1),   mesh(:, 4) - roiBox(2) + 1 + c0(2) - rowNum];
                    
                    % Mesh Rotation with an angle theta
                    xMesh1Rot=zeros(length(xy2Temp(:,1))); % preallocation
                    yMesh1Rot=zeros(length(xy2Temp(:,1)));
                    for meshIdx = 1: length(xy1Temp(:,1))
                        [mesh1Rot] = (R*[xy1Temp(meshIdx,1); xy1Temp(meshIdx,2)])';
                       
                        xMesh1Rot(meshIdx, 1) = mesh1Rot(1, 1);
                        yMesh1Rot(meshIdx, 1) = mesh1Rot(1, 2);
                    end
                   
                    xMesh2Rot=zeros(length(xy2Temp(:,1))); % preallocation
                    yMesh2Rot=zeros(length(xy2Temp(:,1)));
                    for meshIdx = 1: length(xy2Temp(:,1))
                        [mesh2Rot] = (R*[xy2Temp(meshIdx,1); xy2Temp(meshIdx,2)])';
                        xMesh2Rot(meshIdx, 1) = mesh2Rot(1, 1);
                        yMesh2Rot(meshIdx, 1) = mesh2Rot(1, 2);
                    end
                    
                     xy1 = [xMesh1Rot yMesh1Rot];
                     xy2 = [xMesh2Rot yMesh2Rot];
                    
                     % Image rotation
                     imRot = imrotate(imCent, -thetaRot*180/pi); 
                      
                     rotImRow = size(imRot,1);
                     rotImCol = size(imRot,2);
                     
                     xy1ForIm(:, 1) = xy1(:, 1) + rotImCol/2;
                     xy1ForIm(:, 2) = xy1(:, 2) + rotImRow/2;
                     xy2ForIm(:, 1) = xy2(:, 1) + rotImCol/2;
                     xy2ForIm(:, 2) = xy2(:, 2) + rotImRow/2;
                        
                        widthTemp = lWidth;
                        numExtraPix = 6; % look at the profile till a bit outside the mesh % I think actually the pixels outside of the mesh allways have value zero
                        heightTemp = meshLim(ii,2) - meshLim(ii,1) + numExtraPix;
                        
                        initialPointX = rotImCol/2-lWidth/2;
                        initialPointY = rotImRow/2-heightTemp/2;
                        roiBoxTemp =  [initialPointX initialPointY widthTemp heightTemp]; % [xmin ymin width height]

                        % Crop image 
                        imCrop = imcrop(imRot, roiBoxTemp);
                    
                        % Returns pixel intensity values, where n specifies the number of points to include.
                        xHist = size(imCrop,1)/2 - heightTemp/2 :histStep:size(imCrop,1)/2 + heightTemp/2 ;
                        intProf = improfile( imCrop, repmat(size(imCrop,2)/2, 1, length(xHist)), xHist );
                        
%                         % debug: plot intensity profiles
%                         if ii==1 || ii==101 || ii==201
%                             debugFigh=figure;
%                         end
%                         try 
%                             debugFigh.BeingDeleted;
%                         catch
%                             debugFigh=figure;
%                         end
%                         figure(debugFigh)
%                         if ii > 100
%                             if ii > 200
%                                 subplot(10,10,ii-200)
%                             else
%                                 subplot(10,10,ii-100)
%                             end
%                         else
%                             subplot(10,10,ii)
%                         end
%                         plot(xHist,intProf) % or instead run [diamFWHM(ii), ~, ~,~,~, ~, ~]= getFWXM(xHist, intProf', 'searchStep', MINSTEP,'FigZ',debugFigh);
                        
                        slicePosAll{ii}=intProf;
                        MINSTEP=0.01 ;% [pixel]INTERPOLATION STEP - dont expect precision better than MINSTEP*30NM!
        
                        % diamFWHM ------------------------------------------
                        [diamFWHM(ii), ~, ~,~,~, ~, ~]= getFWXM(xHist, intProf', 'searchStep', MINSTEP);
                        
                        % save intensity profiles for later use.
                        intProfCelltmp{1,ii}=xHist;
                        intProfCelltmp{2,ii}=intProf';
                        
                    else
                            intProf = [];
                            slicePosAll{ii} = intProf;
                            if ii == 1
                                diamFWHM(ii) = 0;
                            else
                                diamFWHM(ii) = diamFWHM(ii-1);
                            end
                    end
                end
                
                %% Find maxima and minimum (constriction site) in cell width
                % to find the max diameter - needs to be robust so use spline smoothing
                % Use diamMesh or diamFWHM
                nBreak = 5; % Determines complexity of splinefit to width profile
                ppDiam = splinefit(lengthOut, diamFWHM, nBreak);
                diamSmth = ppval(ppDiam, lengthOut);

                % Find the first local maximum on the left
                [posLeftStart,valLeftStart,idxLeftStart]   =  getLeft(lengthOut ,diamSmth);
                [posRightEnd,valRightEnd,idxRightEnd] = getRight(lengthOut ,diamSmth);
                if isempty(posLeftStart)|| isempty(posRightEnd)
                  maxD=[];
                  minD=[];
                else
                  % then on the right.
                  % this marks the limits of the search area
                  % find the absolute minimum
                  [minVal, idxMin ] = min(diamFWHM(idxLeftStart:idxRightEnd));
                  idxMin = idxMin+idxLeftStart-1;
                  minD = [lengthOut(idxMin), minVal];
                  %find the absolute maximum to the the left and right of the minima
                  %left
                  [maxVal, idxMax] = max(diamSmth(idxLeftStart:idxMin));
                  idxMax= idxMax + idxLeftStart-1;
                  maxD(1,:) = [lengthOut(idxMax), maxVal];
                  %Right
                  [maxVal, idxMax] = max(diamSmth(idxMin:idxRightEnd));
                  idxMax= idxMax +idxMin -1;
                  maxD(2,:) = [lengthOut(idxMax), maxVal];
                end
                
                % Grab the central position where there is the min diam
                c0 = c0AllPos(idxMin, :);
                
                %% Save measurements in cellStruct
                tmpIntProf=intProfCelltmp{:,idxMin};
                cellStruct.MTSprofile=tmpIntProf; % save intensity profile at "constriction site"
             
                cellStruct.bactLength = bactLength;
                cellStruct.origCropImg = roiOrigIm;
                cellStruct.cellIdx = reg;

                cellStruct.polarity = 0;
                cellStruct.box = roiBox;
                cellStruct.ancestors = [];
                cellStruct.descendants = [];
                cellStruct.divisions = [];
                cellStruct.stage = [];
                cellStruct.contour = cCell;
                cellStruct.poleIdx = [pole1PosIdx pole2PosIdx];
                cellStruct.polarity = 0; % algorithm 1 cannot divide cells unless this is interpoutline
                cellStruct.centroid = cellCentroid;
                cellStruct.birthframe = frame;
                if isnan(diamFWHM(1)) % check for nans, but only near the ends of the cell
                    diamFWHM(1)=0;
                    if isnan(diamFWHM(2)) % would there be a more general approach?
                        diamFWHM(2)=0;
                    end
                end
                if isnan(diamFWHM(end))
                    diamFWHM(end)=0;
                    if isnan(diamFWHM(end-1))
                        diamFWHM(end-1)=0;
                    end
                end
                cellStruct.diamFWHM=[lengthOut'; diamFWHM]; % Save the diameter profile over the cell length measured by FWHM. For SA and Volume calculation.
                cellStruct = getSurfAndVolFromCellInfo(cellStruct);%add length area surface volume
                cellStruct.cellAssociatedNum = cellIdx;
                
                cellStruct.maxD = min(maxD(:, 2));
                cellStruct.minD = minD(1, 2);
                cellStruct.waistD = minD(1, 2)/min(maxD(:, 2));
                cellStruct.c0 = c0;
                cellStruct.mesh = mesh;
                

                cell0 = cell0+1;

                if cellIdx == 0
                    cellList{frame, cell0} = cellStruct;
                else
                    cellList = cellStruct;
                end

                display(['fitting region ' num2str(reg) ' - passed and saved as cell0 ' num2str(cell0)])
            else
                % if the cell0 is not on the image border OR it did not pass the
                % quality test - split it
                reason = 'unknown';
                if p.getmesh && length(mesh)<=4, reason = 'bad region shape quality'; end
                if min(cCell(:,1))<=p.imRoi(1), reason = 'cell0 on x=0 boundary'; end
                if min(cCell(:,2))<=p.imRoi(2), reason = 'cell0 on y=0 boundary'; end
                if max(cCell(:,1))>=min(imsizes(1,2), (p.imRoi(1)+p.imRoi(3)-1)), reason = 'cell0 on x=max boundary'; end
                if max(cCell(:,2))>=min(imsizes(1,1), (p.imRoi(2)+p.imRoi(4)-1)), reason = 'cell0 on y=max boundary'; end
                display(['fitting region ' num2str(reg) ' - quality check failed - ' reason])
            end
          %% Plot measured widths vs position along the centreline
          figure(getMeshFig),
            hold off,
            subplot(1,2,2)
            plot(lengthOut,diamMesh, 'k');
            hold all;     
            plot(lengthOut,diamFWHM, 'r');
            plot(lengthOut,diamSmth, 'b');
            plot(posLeftStart,valLeftStart,'mo');
            plot(posRightEnd,valRightEnd,'mo');
            plot(minD(:,1),minD(:,2),'ro');
            plot(maxD(:,1),maxD(:,2),'ko');
            xlabel('Length [pixel]')
            ylabel('Diameter [pixel]')
            legend('mesh not smooth', 'FWHM not smooth', 'FWHM smooth','max left pos', 'max right pos', 'min pos',  'max pos')
            legend off
            drawnow
            %% Save getMeshFig
            mkdir([PathNameResults '\contour\bact' num2str(cellIdx) ])
            if cellIdx == 0
                saveas(gcf, [PathNameResults '\contour\bact' num2str(cellIdx) '\bact' num2str(reg) 'frame' num2str(frIdx) '.fig']);
                print(gcf, '-dtiffn',[PathNameResults '\contour\bact' num2str(cellIdx)  '\bact' num2str(reg) 'frame' num2str(frIdx) '.tif'] );
            else
                saveas(gcf, [PathNameResults '\contour\bact' num2str(cellIdx) '\bact' num2str(cellIdx) 'frame' num2str(frIdx) '.fig']);
                print(gcf, '-dtiffn',[PathNameResults '\contour\bact' num2str(cellIdx)  '\bact' num2str(cellIdx) 'frame' num2str(frIdx) '.tif'] );
            end
        end
    end
end
%------------------------------------------------------
function [se, maskdx, maskdy] = initP()
  % initModel
  
  se = strel('arb',[0 1 0;1 1 1;0 1 0]); % erosion mask, can be 4 or 8 neighbors
  maskdx = fliplr(fspecial('sobel')'); %[-1 0 1; -2 0 2; -1 0 1]/2; % masks for computing x & y derivatives
  maskdy = fspecial('sobel');%[1 2 1; 0 0 0; -1 -2 -1]/2;
end


%% Cell statistics and additional data
function [meshLim] = getCellDist(mesh,c0,dY,dX)

%Get distance above
%rotate the gradient 90degerees to the cLine
perpLine_dX = -1*dY;
perpLine_dY = dX;
perpLine_m = perpLine_dY/perpLine_dX;
meshAbove= getDistLinePoly(mesh(:,1:2),c0,perpLine_m);
%Get distance below
meshBelow = getDistLinePoly(mesh(:,3:4),c0,perpLine_m);
meshLim = [-meshBelow, meshAbove];
end
%---------------------------------------------------
function distMax = getDistLinePoly(polyLine,c0,m)
%%DEBUG
%global cLineSmooth;
%%DEBUG
CELLRADMAX = 2000;
%Formula for (x1,y1) distance L along a straight line,
% given (x0,y0) and m 
y_along = @(y0,L,m,dirn) y0 +(dirn/abs(dirn))*m*L/sqrt(1+m^2);
x_along = @(x0,L,m,dirn) x0 +(dirn/abs(dirn))*L/sqrt(1+m^2);

%find intersection of line centred at c0(2,:) perpendicular to the centreline
%and the top mesh
%find the intersection
foundIsec=false;
distMax = [];
lineEnd = CELLRADMAX;
nMax = 2;
ii=1;
while foundIsec ==false
   X_line(1,:) = [x_along(c0(1),lineEnd,m,-1), y_along(c0(2),lineEnd,m,-1)];
   X_line(2,:) = [x_along(c0(1),lineEnd,m,+1), y_along(c0(2),lineEnd,m,+1)];
   [xIsec,yIsec] = intersections(polyLine(:,1),polyLine(:,2),X_line(:,1),X_line(:,2));
   %%DEBUG
   %hold off;
   %plot(cLineSmooth(:,1),cLineSmooth(:,2));
   %hold all;
   %plot(X_line(:,1),X_line(:,2));
   %plot(polyLine(:,1),polyLine(:,2));
   %if ~isempty(xIsec)
   %  plot(xIsec,yIsec,'ko');
   %end
   %axis equal
   %pause
   %%DEBUG

   if ~isempty(xIsec)
     distMax =sqrt((xIsec-c0(1)).^2+(yIsec-c0(2)).^2);
     distMax = min(distMax);
     foundIsec=true;
   else%increase the length of the search line (max cell radius) if no intersection found
     lineEnd = lineEnd*2;
     if ii ==nMax
       %Intersection between mesh and perpendicular to centreline not found!
       distMax=NaN;
       foundIsec=true;
     else
       ii=ii+1;
     end
   end
end
end

%---------------------------------------------------------
function [posLeftMax,valLeftMax,idxLeftMax] = getLeft(x,y)
nX = numel(x);
ii =2;
foundLeftPk = false;
idxLeftMax=[];
valLeftMax=[];
posLeftMax=[];
while ~foundLeftPk && ii<=nX-1
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<0)%absolute maximum
    idxLeftMax = ii;
    valLeftMax = y(ii);
    posLeftMax=x(ii);
    foundLeftPk=true;
 % elseif (y(ii)-y(ii-1)<0)%if we're descending we're past the peak
 %   idxLeftMax = 1;
 %   valLeftMax = y(1);
 %   posLeftMax=x(1);
 %   foundLeftPk=true;
  else
    ii=ii+1;
  end
end
end

%---------------------------------------------------------
function [posRightMax,valRightMax,idxRightMax] = getRight(x,y)
nX = numel(x);
ii =nX-1;
foundRightPk = false;
idxRightMax=[];
valRightMax=[];
posRightMax=[];
while ~foundRightPk && ii>=2
  if (y(ii)-y(ii-1)>0)&&(y(ii+1)-y(ii)<0)
    idxRightMax = ii;
    valRightMax = y(ii);
    posRightMax=x(ii);
    foundRightPk=true;
%  elseif (y(ii)-y(ii-1)>0)%if we're descending we're past the peak
%    idxRightMax =nX;
%    valRightMax = y(nX);
%    posRightMax=x(nX);
%    foundRightPk=true;
  else
    ii=ii-1;
  end
end


end

%-------------------------------------------------------------
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