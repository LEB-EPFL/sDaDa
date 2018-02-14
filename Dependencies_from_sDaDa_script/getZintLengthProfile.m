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

function [cellInfo, ZintAll]=getZintLengthProfile(cellInfo)
% Function to extract Z intensity profile from saved figure in case
% it wasn't saved in cellInfo.

%% initialize
PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160621 DC FtsZ ML159 35min\';
FOV_base='series';
FOVs={'7\Results\160621 DC FtsZ ML159 35min' '8\Results\160621 DC FtsZ ML159 35min' '10\Results\160621 DC FtsZ ML159 35min' '11\Results\160621 DC FtsZ ML159 35min' '13\Results\160621 DC FtsZ ML159 35min' '15\Results\160621 DC FtsZ ML159 35min'};
% % name correction Ambroise
% for fovIdx=1:length(FOVs)
%     FOVs{1,fovIdx}=[FOVs{1,fovIdx} '\Results\160613 ML2159serie' FOVs{1,fovIdx}];
% end
FOVcells=[8 5 6 9 12 8];
diamFWHMAll{sum(FOVcells),75}=[];
numframes=75;
%% loop
cellIdx2=1;
for fovIdx=1:length(FOVs)
    numCells=FOVcells(fovIdx);
    for cellIdx1=1:numCells
        for frIdx=1:numframes
            % set filename
            figFile=[PathNameResults, FOV_base, FOVs{1,fovIdx}, '\contour\bact', num2str(cellIdx1), '\Z_bact', num2str(cellIdx1), 'frame', num2str(frIdx), '.fig'];
            % open figure
            if exist(figFile,'file')
            uiopen(figFile,1)
            % extract data
            fig=gcf;
            Zint=extrZintFromFigure(fig);
            close(fig)
            % save data
            cellInfo{frIdx,1}{1,cellIdx2}.ZintLength=Zint;
            ZintAll{cellIdx2,frIdx}=Zint;
            end
        end
        cellIdx2=cellIdx2+1;
    end
end
end

function Zint=extrZintFromFigure(fig)
% get data from the figure
axesObjs = get(fig, 'Children');
dataObjs = get(axesObjs(2), 'Children');
try
    xdata = get(dataObjs(5), 'XData');
    ydata = get(dataObjs(5), 'YData');
catch % this can happen when the text is not on the plot.
    xdata = get(dataObjs(4), 'XData');
    ydata = get(dataObjs(4), 'YData');
end
Zint=[xdata;ydata];
end
