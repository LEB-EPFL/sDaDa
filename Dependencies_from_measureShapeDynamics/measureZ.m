function cellList = measureZ(imgZ, cellList, ~, roiCor, cellIdx, frIdx, PathNameResults,FtsZFig)
% cellList = measureZ(imgZ, cellList, ~, roiCor, cellIdx, frIdx, PathNameResults,FtsZFig)
% INPUT:
% imgZ: image of the cell's ROI in the divisome channel
% cellList: cellInfo structure for the current frame
% roiCor: Drift correction for ROI
% cellIdx: cell number, used for plotting and saving figures (naming)
% frIdx:  frame number, used for plotting and saving figures (naming)
% PathNameResults: path for saving figures
% FtsZFig: figure handle for the figure where the measured diameter, the
% contour, mesh etc... are plotted.
% OUTPUT:
% cellList:
% A figure with the results will be plotted and saved
% This function is based on getMeshAlgForSIM.
% Its purpose is to evaluate/measure the divisome intensity images for dual
% color SIM analysis. It will use the information on cell location and
% shape from the MTS channel to measure the integrated intensity profile
% along the length of the cell. Where the intensity is maximal (at the
% divisome), a cross-section of the intensity will be taken and its with
% measured. When there is a peak around midcell that is significantly
% bigger than all other peaks, the divisome will be considered as
% "assembled", reflected in the cellList.Zassembled value.
% 
% Author: Aster Vanhecke

%% Grab input and intialize parameters
mesh=cellList.mesh;
roiOrigIm=imgZ;
histStep = 0.5;
pixSize=30;
lStep = 2; % change
lWidth = 5;% change
cLine=getCLine(cellList);

%% Crop image, change coordinates of mesh and cLine
roiBox=cellList.box+[-5 -5 10 10];
roiBox(1)=roiBox(1)+roiCor(1);
roiBox(2)=roiBox(2)+roiCor(2);
roiOrigIm=imcrop(roiOrigIm,roiBox);
mesh(:,1)=mesh(:,1)-roiBox(1)+roiCor(1);
mesh(:,2)=mesh(:,2)-roiBox(2)+roiCor(2);
mesh(:,3)=mesh(:,3)-roiBox(1)+roiCor(1);
mesh(:,4)=mesh(:,4)-roiBox(2)+roiCor(2);
cLine(:,1)=cLine(:,1)-roiBox(1)+roiCor(1);
cLine(:,2)=cLine(:,2)-roiBox(2)+roiCor(2);

%% Plot subplot with contour and cLine on the divisome channel image
figure(FtsZFig)
subplot(2,2,1),
hold off,
imshow(roiOrigIm, 'InitialMagnification',500, 'Colormap', jet, 'DisplayRange', [0 max(max(roiOrigIm))+1e-10]),
hold on,
plot(mesh(:,1),mesh(:,2)),
plot(mesh(:,3),mesh(:,4)),
plot(cLine(:,1),cLine(:,2)),
title(['Cell ' num2str(cellIdx) ', frame ' num2str(frIdx)])

%% Step along centerline and measure intensity
% Prepare length steps
[~, seglen]= arclength(cLine(:,1), cLine(:,2), 'lin');
len=[0; cumsum(seglen)];

l0 = 0:lStep:max(len);
nStep = numel(l0)-1;

% Step and Measure
intensities=zeros(nStep,1);
for ii = 1:nStep
    [intensities(ii), ~]=rotateNInt(ii);
end

%% Find intensity peak, plot intensity along centerline
figure(FtsZFig),
subplot(2,2,3:4),
hold off,
plot(l0(1:end-1).*pixSize,intensities), title(['Z intensity projection over cell length for cell ' num2str(cellIdx) ' in frame ' num2str(frIdx)])
maxZ=find(intensities==max(intensities));
hold on,
thresh=graythresh(intensities);
plot([0 l0(end-1)].*pixSize,[thresh thresh],'g')
plot([0 l0(end-1)].*pixSize,[intensities(maxZ)/3 intensities(maxZ)/3], 'r')

%% define midcell region.
midCellStart=0.3*(nStep+1);
midCellStop=0.7*(nStep+1);
plot([l0(round(midCellStart)) l0(round(midCellStart)) l0(round(midCellStop)) l0(round(midCellStop))].*pixSize,[0 1 1 0],'b')
drawnow
cellList.ZintLength=[l0(1:end-1).*pixSize;intensities'];

%% Check if divisome is assembled
Zassembled=false;
if midCellStart<maxZ && maxZ < midCellStop % If the highest peak is at midcell
    % And that the intensities outside the midcell region do not exceed one third of the highest (midcell) peak
    if ~sum(intensities(1:round(midCellStart))>intensities(maxZ)/3) && ~sum(intensities(round(midCellStop):end)>intensities(maxZ)/3)
        Zassembled=true;
    end
end
cellList.Zassembled=Zassembled;

if Zassembled % Write result on fig
    text(10, 0.8, 'The Z-ring has assembled!')
else
    text(10, 0.8, 'The Z-ring has NOT assembled!')
end

% Record intensity fraction at midcell region
midcellInt= sum( intensities(round(midCellStart):round(midCellStop)) )/sum(intensities);
cellList.midcellInt=midcellInt;

%% Measure FWHM at intensity peak, plot intensity cross-section
[~, xHist, intProf]=rotateNInt(maxZ);
cellList.diamZ=[];
figure(FtsZFig),
subplot(2,2,2)
[diamZ, ~, ~,~,~]=getFWXM(xHist, intProf','FigZ',FtsZFig);
cellList.diamZ=diamZ*pixSize; % save diameter of Zring
cellList.Zpos=l0(maxZ)/len(end); % save position of Zring
cellList.Zprofile=[xHist' intProf];

% save figure
saveas(FtsZFig, [PathNameResults '\contour\bact' num2str(cellIdx) '\Z_bact' num2str(cellIdx) 'frame' num2str(frIdx) '.fig']);
print(FtsZFig, '-dtiffn',[PathNameResults '\contour\bact' num2str(cellIdx)  '\Z_bact' num2str(cellIdx) 'frame' num2str(frIdx) '.tif'] );

%% Function for image rotation and taking intensity profile
function [Int, xHist, intProf]=rotateNInt(ii)
    % Find current mesh centre position, and +- lWidth/2
    % Probably best to smooth this!
    cPt(:,1) = interp1(len,cLine(:,1),[l0(ii), l0(ii)+lWidth/2, l0(ii)+lWidth]');
    cPt(:,2) = interp1(len,cLine(:,2),[l0(ii), l0(ii)+lWidth/2, l0(ii)+lWidth]');
    c0Temp = cPt(2,:);%this is the midpoint
    
    % Get angle
    dX = cPt(3,1)-cPt(1,1);
    dY = cPt(3,2)-cPt(1,2);
    theta = atan2(dY,dX);
    
    % Get meshMin Max distances
    % meshLim = [-meshBelow, meshAbove];
    meshLim(ii,:) = getCellDist(mesh,c0Temp,dY,dX);
    
    if ~any(isnan(meshLim(ii, :)))
        rowNum = size(roiOrigIm, 1);
        colNum = size(roiOrigIm, 2);
        
        c0(1)= c0Temp(1);
        c0(2)= rowNum - c0Temp(2);
        
        distYPointFirstRow = rowNum - c0(2);
        distXPointLastCol = colNum - c0(1);
        
        % Image translation --> c0 should be the center of the image,
        % for the rotation
        % NOTE I changed floor to round my function trans
        if distYPointFirstRow*2 >= rowNum && c0(1)*2 >= colNum
            imCent = zeros(round(distYPointFirstRow*2)+1, round(c0(1)*2)-1);
            imCent(1: rowNum, 1:colNum) = roiOrigIm;
        elseif distYPointFirstRow*2 < rowNum && c0(1)*2 < colNum
            if round( c0(2) )*2 - 1 > rowNum % this one was used
                imCent = zeros(round(c0(2)*2)-1, round(distXPointLastCol*2)+1);
                newColNum = round(distXPointLastCol*2)+1;
                newRowNum = round( c0(2) )*2 - 1;
                imCent(newRowNum - rowNum + 1 : newRowNum,newColNum - colNum + 1 : newColNum) = roiOrigIm;
            else
                imCent = zeros(round(c0(2)*2)+1, round(distXPointLastCol*2)+1);
                newColNum = round(distXPointLastCol*2)+1;
                newRowNum = round( c0(2) )*2 +1;
                imCent(newRowNum - rowNum + 1 : newRowNum , newColNum - colNum + 1 : newColNum) = roiOrigIm;
            end
        elseif distYPointFirstRow*2 < rowNum && c0(1)*2 >= colNum
            if round( c0(2) )*2 - 1 > rowNum
                imCent = zeros(round(c0(2)*2)-1, round(c0(1)*2)-1);
                newRowNum = round(c0(2)*2)-1;
                imCent(newRowNum - rowNum + 1 : newRowNum, 1:colNum) = roiOrigIm;
            else
                imCent = zeros(round(c0(2)*2)+1, round(c0(1)*2)-1);
                newRowNum = round(c0(2)*2)+1;
                imCent(newRowNum - rowNum + 1 : newRowNum, 1:colNum) = roiOrigIm;
            end
        elseif distYPointFirstRow*2 >= rowNum && c0(1)*2 < colNum
            imCent = zeros(round(distYPointFirstRow*2)+1, round(distXPointLastCol*2)+1);
            newColNum = round(distXPointLastCol*2)+1;
            imCent(1: rowNum, newColNum - colNum + 1 : newColNum) = roiOrigIm;
        end
        
        %% Rotate parallel to xAxis
        thetaRot = -theta;
        R = [cos(thetaRot),-sin(thetaRot );...
            sin(thetaRot ),cos(thetaRot)];
        
        %% Mesh Translation to center in c0
        xy1Temp = [mesh(:, 1) + 1 - c0(1),   mesh(:, 2) + 1 + c0(2) - rowNum];
        xy2Temp = [mesh(:, 3) + 1 - c0(1),   mesh(:, 4) + 1 + c0(2) - rowNum];
        
        %% Mesh Rotation by angle theta
        xMesh1Rot=zeros(size(xy1Temp(:,1),1)); %preallocation
        yMesh1Rot=zeros(size(xy1Temp(:,1),1)); %preallocation
        
        for meshIdx = 1: length(xy1Temp(:,1))
            [mesh1Rot] = (R*[xy1Temp(meshIdx,1); xy1Temp(meshIdx,2)])';
            
            xMesh1Rot(meshIdx, 1) = mesh1Rot(1, 1);
            yMesh1Rot(meshIdx, 1) = mesh1Rot(1, 2);
        end
        
        xMesh2Rot=zeros(size(xy2Temp(:,1),1)); %preallocation
        yMesh2Rot=zeros(size(xy2Temp(:,1),1)); %preallocation
        for meshIdx = 1: length(xy2Temp(:,1))
            [mesh2Rot] = (R*[xy2Temp(meshIdx,1); xy2Temp(meshIdx,2)])';
            xMesh2Rot(meshIdx, 1) = mesh2Rot(1, 1);
            yMesh2Rot(meshIdx, 1) = mesh2Rot(1, 2);
        end
        
        %% Image rotation by angle theta
        imRot = imrotate(imCent, -thetaRot*180/pi);
        
        rotImRow = size(imRot,1);
        rotImCol = size(imRot,2);
        
        numExtraPix = 6; % look at the profile till a bit outside the mesh
        heightTemp = meshLim(ii,2) - meshLim(ii,1) + numExtraPix;
        
        initialPointX = rotImCol/2-lWidth/2;
        initialPointY = rotImRow/2-heightTemp/2;
        roiBoxTemp =  [initialPointX initialPointY lWidth heightTemp]; % [xmin ymin width height]
        
        % Crop image
        imCrop = imcrop(imRot, roiBoxTemp);
        
        %% Returns pixel intensity values, where n specifies the number of points to include.
        xHist = size(imCrop,1)/2 - heightTemp/2 :histStep:size(imCrop,1)/2 + heightTemp/2 ;
        intProf = improfile( imCrop, repmat(size(imCrop,2)/2, 1, length(xHist)), xHist );
        intProf(isnan(intProf))=0;
        Int=sum(intProf);
        
    else
        Int=0;
        intProf=0;
        xHist=0;
    end
end
end

function [meshLim] = getCellDist(mesh,c0,dY,dX)

%Get distance above
%rotate the gradient 90degerees to the cLine
perpLine_dX = -dY;
perpLine_dY = dX;
perpLine_m = perpLine_dY/perpLine_dX;
meshAbove= getDistLinePoly(mesh(:,1:2),c0,perpLine_m);
%Get distance below
meshBelow = getDistLinePoly(mesh(:,3:4),c0,perpLine_m);
meshLim = [-meshBelow, meshAbove];
end

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



