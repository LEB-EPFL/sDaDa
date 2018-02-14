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
function [cellInfo, diamFWHMAll]=getDiamFWHMlengthProfile(cellInfo)
% Function to extract diameter (by FWHM) profile from saved figure in case
% it wasn't saved in cellInfo.

%% initialize
PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160613 DC ML2159 t0 80min\Results\';
FOV_base='serie';
FOVs={'2\Results\160613 ML2159serie2' '4\Results\160613 ML2159serie4' '5\Results\160613 ML2159serie5' '6\Results\160609 ML2159serie6' '7\Results\160613 ML2159DC t80'};
% % name correction Ambroise
% for fovIdx=1:length(FOVs)
%     FOVs{1,fovIdx}=[FOVs{1,fovIdx} '\Results\160613 ML2159serie' FOVs{1,fovIdx}];
% end
FOVcells=[14 12 10 13 16];
diamFWHMAll{sum(FOVcells),75}=[];
numframes=75;
%% loop
cellIdx2=1;
for fovIdx=1:length(FOVs)
    numCells=FOVcells(fovIdx);
    for cellIdx1=1:numCells
        for frIdx=1:numframes
            % set filename
            figFile=[PathNameResults, FOV_base, FOVs{1,fovIdx}, '\contour\bact', num2str(cellIdx1), '\bact', num2str(cellIdx1), 'frame', num2str(frIdx), '.fig'];
            % open figure
            if exist(figFile,'file')
            uiopen(figFile,1)
            % extract data
            fig=gcf;
            diamFWHM=extrDiamFWHMFromFigure(fig);
            close(fig)
            try
                if max(diamFWHM(2,:)) < 2
                    diamFWHM=[];
                end
                % replace nans near the end of the cell by 0
                if isnan(diamFWHM(2,1)) || diamFWHM(2,1)==0 % check for nans, but only near the ends of the cell
                    diamFWHM(2,1)=0;
                    if isnan(diamFWHM(2,2)) % would there be a more general approach?
                        diamFWHM(2,2)=0;
                    end
                end
                if isnan(diamFWHM(2,end)) || diamFWHM(2,end)==0
                    diamFWHM(2,end)=0;
                    if isnan(diamFWHM(2,end-1))
                        diamFWHM(2,end-1)=0;
                    end
                end
            catch % if diamFWHM=[] this probably didnt work
                size(diamFWHM)
                warning('line 50:')
            end
            % save data
            cellInfo{frIdx,1}{1,cellIdx2}.diamFWHM=diamFWHM;
            diamFWHMAll{cellIdx2,frIdx}=diamFWHM;
            end
        end
        cellIdx2=cellIdx2+1;
    end
end
end

function diamFWHM=extrDiamFWHMFromFigure(fig)
% get data from the figure
axesObjs = get(fig, 'Children');
try
    dataObjs = get(axesObjs(1), 'Children');
    try
        xdata = get(dataObjs(6), 'XData');
        ydata = get(dataObjs(6), 'YData');
    catch
        try
            dataObjs = get(axesObjs(2), 'Children');
            xdata = get(dataObjs(6), 'XData');
            ydata = get(dataObjs(6), 'YData');
        catch
            dataObjs = get(axesObjs(3), 'Children');
            xdata = get(dataObjs(6), 'XData');
            ydata = get(dataObjs(6), 'YData');
        end
    end
    diamFWHM=[xdata;ydata];
catch
    axesObjs = get(fig, 'Children')
    diamFWHM=[];
end
end
