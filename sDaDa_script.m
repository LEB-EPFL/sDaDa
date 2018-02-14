% This is a scrip for analysing cell images
% - segmentation
% - cell parameter measured (diameter, length, pole positions, etc)
%
% The main functions in the script returns in OUTPUT:
% - imageStack: grayscaled, normalized imageStack (for debugging/testing)
% - cellInfo: a cell cointaining all the info af all the bacteria at each
%             frame
% - twMatrix: the matrix with the time and the waist diameter for each
%             bacteria and frame (see function pltSIM_diameter.m)
%                   %       cell1  cell2 ... celln
%                     fr1   t1 w1  t2 w2 ... tn wn
%                     fr2
% - tlMatrix: the matrix with the time and the length for each
%             bacteria and frame (see function pltSIM_diameter.m)
%                   %       cell1  cell2 ... celln
%                     fr1   t1 l1  t2 l2 ... tn ln
%                     fr2
% - Tvar: Tc, Tg from the fitting from both feingold and 2D model (see
%         function pltSIM_diameter.m)  Tvar = [Tc' Tg' TcF' TgF'];
% - Lvar: Lc, Lg from the fitting from both feingold and 2D mode (see
%         function pltSIM_diameter.m)  Lvar = [Lc' Lg' LcF' LgF'];

%%
tic
%% PARAMETERs
tInt = 5; % Time between consecutive frames [min]
T0 = 30; % Starting time with respect to synchrony [min]
pixSize = 30; % Pixel size [nm]
% Use this line to manually set parameters and resultsfolder for rerunning pltSIM_diameter
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

%% Do the segmentation of the SIM data and measure some quantities(waist diameter, contour, etc)
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

%% Plot diameter vs time, length vs time etc and do the fit
diamThresh=90; smthThresh=0.92;
[twMatrix, tlMatrix, Tvar, Lvar, BIC, p1C0, meanR, twZMat, tZMat, fitResultF, fitResultX, fitResultFz, fitResultXz, tDMatrix] = analyzeShapeDynamics(cellInfo, tInt, T0, pixSize, divededVarEachFrm,PathNameResults, diamThresh,smthThresh);
ElMat=getElMat(tlMatrix);

%% Save variable (change this)
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