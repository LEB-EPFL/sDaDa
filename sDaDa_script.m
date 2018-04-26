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

% This is a scrip for analysing cell images
% - segmentation
% - cell parameter measured (diameter, length, pole positions, etc)
%
% The main functions in the script returns in OUTPUT:
% - imageStack: grayscaled, normalized imageStack (for debugging/testing)
% - cellInfo: a cell cointaining all the info af all the bacteria at each
%             frame
% - twMatrix: the matrix with the time and the waist diameter for each
%             bacteria and frame (see function analyzeShapeDynamics.m)
%                   %       cell1  cell2 ... celln
%                     fr1   t1 w1  t2 w2 ... tn wn
%                     fr2
% - tlMatrix: the matrix with the time and the length for each
%             bacteria and frame (see function analyzeShapeDynamics.m)
%                   %       cell1  cell2 ... celln
%                     fr1   t1 l1  t2 l2 ... tn ln
%                     fr2
% - Tvar: Tc, Tg from the fitting from both feingold and 2D model (see
%         function analyzeShapeDynamics.m)  Tvar = [Tc' Tg' TcF' TgF'];
% - Lvar: Lc, Lg from the fitting from both feingold and 2D mode (see
%         function analyzeShapeDynamics.m)  Lvar = [Lc' Lg' LcF' LgF'];

%%
tic
printLicense_sDaDa
%% PARAMETERs
tInt = 5; % Time between consecutive frames [min]
T0 = 30; % Starting time with respect to synchrony [min]
pixSize = 30; % Pixel size [nm]
% Use this line to manually set parameters and resultsfolder for rerunning analyzeShapeDynamics
%tInt=5;T0=22;pixSize=30;PathNameResults='\\Data\SIM data\Analysis\170307_ML2159_label_10\Results\ExperimentName'; 

%% File name and path
% Naming convention
% membrane channel: PathName\dataSet_FOVnr_MTS.tif
% divisome channel: PathName\dataSet_FOVnr_FtsZ.tif % FtsZ can be changed
% if the definition of FileNameZ is changed accordingly (line 48)
% Results are saved in PathName\Results\dataSetFOVnr\
dataSet='160713 DC WT t0 30_'; % change this
FOVnr='03'; % change this
FileName = [dataSet FOVnr '_MTS.tif'];
PathName = '.\Example dataset\'; % change this, include final "\"
FileBaseName = [dataSet FOVnr ]; % names results folder
FileNameZ = [dataSet FOVnr '_FtsZ.tif' ]; % FtsZ channel, comment out for single color analysis
PathNameResults = [PathName 'Results\' FileBaseName];

%% measureShapeDynamics
% Do the segmentation of the SIM data and measure some quantities (waist diameter, contour, etc)
if ~exist('FileNameZ','var')
    [imageStack,  num_frames, cellInfo, divededVarEachFrm] = measureShapeDynamics(FileName, PathName, PathNameResults); % Single color
else
    [imageStack,  num_frames, cellInfo, divededVarEachFrm] = measureShapeDynamics(FileName, PathName, PathNameResults, FileNameZ); %Dual color
end

%% Correct cellInfo and divededVarEachFrm for cells that were incorrectly categorized as not divided.
[cellInfo,divededVarEachFrm] = tooLateCorrect(cellInfo, divededVarEachFrm);

%% Save variables
save([PathNameResults '\cellInfoWithC0.mat'], 'cellInfo')
save([PathNameResults '\cellDivVarWithC0.mat'], 'divededVarEachFrm')

%% analyzeShapeDynamics
% Plot diameter vs time, length vs time, perform fits, determine
% constriction onset time, etc.
diamThresh=90; smthThresh=0.92;
[twMatrix, tlMatrix, Tvar, Lvar, BIC, p1C0, meanR, twZMat, tZMat, fitResultF, fitResultX, fitResultFz, fitResultXz, tDMatrix] = analyzeShapeDynamics(cellInfo, tInt, T0, pixSize, divededVarEachFrm,PathNameResults, diamThresh,smthThresh);
ElMat=getElMat(tlMatrix);

%% Save variables
save([PathNameResults '\twMatrix.mat'], 'twMatrix');
save([PathNameResults '\tlMatrix.mat'], 'tlMatrix');
save([PathNameResults '\Tvar.mat'], 'Tvar');
save([PathNameResults '\Lvar.mat'], 'Lvar');
save([PathNameResults '\p1C0.mat'], 'p1C0');
save([PathNameResults '\meanR.mat'], 'meanR');
save([PathNameResults '\twZMat.mat'], 'twZMat');
save([PathNameResults '\tZMat.mat'], 'tZMat');
save([PathNameResults '\LZvar.mat'], 'LZvar');
save([PathNameResults '\fitResultF.mat'], 'fitResultF');
save([PathNameResults '\fitResultX.mat'], 'fitResultX');
save([PathNameResults '\fitResultFz.mat'], 'fitResultFz');
save([PathNameResults '\fitResultXz.mat'], 'fitResultXz');
save([PathNameResults '\ElMat.mat'], 'ElMat');

%% Extract and save Dmax for all cells
DmaxAll=getDmaxAll(cellInfo);
save([PathNameResults '\DmaxAll.mat'],'DmaxAll')
tTemp=T0:tInt:T0+(num_frames-1)*tInt;
Tc=Tvar(:,10);
for cellIdx=size(DmaxAll,1):-1:1
    DmaxAtTcAll(cellIdx,1)=nanmean(DmaxAll(cellIdx,Tc(cellIdx,1)-tInt <= tTemp & tTemp <= Tc(cellIdx,1)+tInt),2);
end
save([PathNameResults '\DmaxAtTcAll.mat'],'DmaxAtTcAll')

%% Extracting Z intensity and width profiles, calculating SAV and saving all of these.
ZintAll=getZintAll(cellInfo);
diamFWHMAll=getDiamFWHMAll(cellInfo);
cellInfo=getSAVforCellInfo(cellInfo,1);
[SAVall, SAall, Vall]=getSAVall(cellInfo);

save([PathNameResults '\ZintAll.mat'], 'ZintAll');
save([PathNameResults '\diamFWHMAll.mat'], 'diamFWHMAll');
save([PathNameResults '\SAVall.mat'], 'SAVall');
save([PathNameResults '\SAall.mat'], 'SAall');
save([PathNameResults '\Vall.mat'], 'Vall');

%%
duration=toc;
fprintf('Finished! Elapsed time is:\n')
breaktime(duration,'disp');