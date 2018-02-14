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

% script to run analyzeShapeDynamics overnight :)
% remember to configure togglePlots
% does:
% load cellInfo and celldivededvar
% set parameters specify! T0 etc
% run analyzeShapeDynamics
% save output
% close figures
% try SAV stuff
% clear and close all
tic
for analysisCycleId=5 % play with this
    %% ------------
    tInt=5;
    pixSize=30;
    switch analysisCycleId
        case 1
            T0=30;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160713 DC WT t0 30min\Results\together_05-12';
        case 2
            T0=22;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160622 DC WT t0 22min\Results\together 3-5_9_10_12-15';
        case 3
            T0=35;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160621 DC FtsZ ML159 35min\together 2_4-11_13_15';
        case 4
            T0=80;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160613 DC ML2159 t0 80min\Results\together 1_2_4-9';
        case 5
            T0=34;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\161005 FtsW Mut 34min\Results\together 1-12\manual Tw 2';
        case 6
            T0=50;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160816 WT DC 1305 FtsWGFP t0 50min\Results\together 01-06 and 10\manual Tw 2';
        case 7
            T0=37;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160926 WT FtsW 18 488 37 min\Results\Together 01-07\manual Tw 2';
        case 8
            T0=36;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160927 DC ML2159 FtsW t0 36min\Results\together_01-10\manual Tw 3';
        case 9 % FOM
            T0=80; tInt=10;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160526 dc fOSFO T080\Results\together 1-9, 11,12,14-17';
        case 10 % FtsN A
            T0=45; tInt=5;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160727 WT FTsN 45min\Results\together 1-11\';
        case 11
            T0=28;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170412 WT DC FtsN t28\Results\together 1-4 11-15\';
        case 12
            T0=50;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160726 ML2159 DC FtsN t0 50min\Results\together 1_6_7_9\';
        case 13 
            T0=25;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160810 ML2159 FtsN t0 25 min\Results\ML2159 FtsN t25 01_06_07_09 together\';
        case 14
            T0=22;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170221 WT PYE T0 22min\Results\together 1-5_7_9_11\Cleaned\';
        case 15 % PYE Mut
            T0=16;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170307 ML2159 DC FtsZ PYE\Results\together 1-3_5_6_8_14_15\';
            
    end
    %%
    load([PathNameResults '\cellInfoWithC0.mat'])
    load([PathNameResults '\cellDivVarWithC0.mat'])
    PathNameResults=[PathNameResults '\rerunForFig'];
    mkdir(PathNameResults);
    cellInfo = rerunGMAFS(cellInfo,PathNameResults);
    save([PathNameResults '\cellInfoWithC0.mat'], 'cellInfo');
    
    %% redo analyzeShapeDynamics since diameter threshold onset might have changed
    [twMatrix, tlMatrix, Tvar, Lvar, BIC, p1C0, meanR, twZMat, tZMat, TZvar, LZvar, fitResultF, fitResultX, fitResultFz, fitResultXz] = analyzeShapeDynamics(cellInfo, tInt, T0, pixSize, divededVarEachFrm,PathNameResults, 90);
    
    %% Save variable (change this)
    save([PathNameResults '\twMatrix.mat'], 'twMatrix');
    save([PathNameResults '\tlMatrix.mat'], 'tlMatrix');
    save([PathNameResults '\Tvar.mat'], 'Tvar');
    save([PathNameResults '\Lvar.mat'], 'Lvar');
%     save([PathNameResults '\p1C0.mat'], 'p1C0');
%     save([PathNameResults '\meanR.mat'], 'meanR');
    save([PathNameResults '\twZMat.mat'], 'twZMat');
    save([PathNameResults '\tZMat.mat'], 'tZMat');
    save([PathNameResults '\LZvar.mat'], 'LZvar');
    save([PathNameResults '\fitResultF.mat'], 'fitResultF');
    save([PathNameResults '\fitResultX.mat'], 'fitResultX');
    save([PathNameResults '\fitResultFz.mat'], 'fitResultFz');
    save([PathNameResults '\fitResultXz.mat'], 'fitResultXz');
    %%
    ElMat=getElMat(tlMatrix);
    save([PathNameResults '\ElMat.mat'], 'ElMat');

    %% Extracting width profiles, calculating SAV and saving all of these.
    try
        diamFWHMAll=getDiamFWHMAll(cellInfo);
        cellInfo=getSAVforCellInfo(cellInfo,1);
        [SAVall, SAall, Vall]=getSAVall(cellInfo);
        save([PathNameResults '\diamFWHMAll.mat'], 'diamFWHMAll');
        save([PathNameResults '\SAVall.mat'], 'SAVall');
        save([PathNameResults '\SAall.mat'], 'SAall');
        save([PathNameResults '\Vall.mat'], 'Vall');
    catch
        disp(['getting and saving SAV failed for dataset ' num2str(analysisCycleId)])
    end
    save([PathNameResults '\cellInfoWithC0.mat'], 'cellInfo');
    clear
    close all
end
%%--------------
duration=toc;
fprintf('Finished! Elapsed time is:\n')
breaktime(duration,'disp');
