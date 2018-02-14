% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Aster Vanhecke 
% 
% Contact:
% e-mail:
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

function cellInfo = rerunGMAFS(cellInfo,PathNameResults)
% function to rerun getMeshAlgForSIM for a dataset
param=initP;
num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
scrsz = get(0,'ScreenSize');
getMeshFig=figure('Position',[scrsz(1,3)/2 50 scrsz(1,3)/2 scrsz(1,4)/2]);
for cellIdx=1:num_cells
    for frIdx=1:num_frames
        skip=false;
        try
            img=cellInfo{frIdx, 1}{1, cellIdx}.origCropImg;
        catch
            warning(['Failed for cell ' num2str(cellIdx) ' in frame ' num2str(frIdx)])
            skip=true;
        end
        if ~skip
            % Expand image using roiBox. Important for coordinates of Mesh etc. 
            box=cellInfo{frIdx, 1}{1, cellIdx}.box;
            imgExpanded=zeros(ceil(box(1:2))+ceil(box(3:4)) + [5 5]);
            imgExpanded( ceil(box(1)):ceil(box(1))+size(img,1)-1 , ceil(box(2)):ceil(box(2))+size(img,2)-1 ) = img;
            cellList = getMeshAlgForSIM(imgExpanded, imgExpanded>0, param, cellIdx, frIdx, PathNameResults,getMeshFig);
            cellInfo{frIdx, 1}{1, cellIdx}.MTSprofile = cellList.MTSprofile;
            cellInfo{frIdx, 1}{1, cellIdx}.diamFWHM = cellList.diamFWHM;
            cellInfo{frIdx, 1}{1, cellIdx}.maxD = cellList.maxD;
            cellInfo{frIdx, 1}{1, cellIdx}.minD = cellList.minD;
            cellInfo{frIdx, 1}{1, cellIdx}.waistD = cellList.waistD;
            cellInfo{frIdx, 1}{1, cellIdx}.c0 = cellList.c0;
        end
    end
end
end

function param=initP()
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
end