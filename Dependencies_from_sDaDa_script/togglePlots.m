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

%% version to only redo Lvar/Tvar
%% To plot or not to plot?
% plotSCWvT=              false; % Single Cell Width vs Time
%     plotSC_fitRes=      false; % plot fit residuals?
% plotSCLvT=              false; % Single Cell Length vs Time
%     plotAllLvT=         false; % All cells Length vs Time
%         plotAllnormLvT= false; % All cells normalized Length vs Time
% plotdWvT=               false; % "constriction": dWidth vs time
% plotBIChist=            false; % histogram of BIC of fits
% plotSCWvTZ=             false; % Single Cell Divisome and inner membrane diameter vs Time 
% plotCurv=               false; % Curvature fit (debug) (slow!)
% plotCurvAll=            false; % All cells Curvature vs normalized time
% plotCentPos=            false; % All cells normalized position of constriction site
% 
% % To fit or not to fit?
% fitSC_Feingold=         true;  % Feingold model to single cell data
% fitSC_Xiaofree=         false; % Xiao model to single cell data
% fitSC_2D=               false; % 2D model to single cell data
% fitSCZ_Feingold=        false; % Feingold model to single cell divisome diameter data
% fitSCZ_Xiaofree=        false; % Xiao model to single cell divisome diameter data
% fitSCmts_Feingold=      false; % Feingold model to single cell diameter data
% fitSCmts_Xiaofree=      false; % Xiao model to single cell diameter data
% fitSC_Xiaofree_Tz=      true;  % Xiao model with fixed Tc to single cell data
% fitSC_Xiaofree_Tc_Tg =  false; % Xiao model with fixed Tc and Tg to single cell data % Overwrites fitSC_Xiaofree_Tz
% fitSC_Feingold_Tz=      false; % Feingold model with fixed Tc to single cell data
% fitSC_XiaoAlpha =       false; % Xiao model with fixed Tc to single cell data
%
% % To calculate or not to calculate?
% calcCurv=               false; % Radius of curvature of cell center line
% calcMTSmid=             false; % Position of central diameter local minimum (on long axis)
% 
% % ensure dependencies
% if ~calcMTSmid
%     plotCurv=           false;
%     plotCurvAll=        false;
% end
% if ~calcCurv
%     plotCurv=           false;
%     plotCurvAll=        false;
% end

%%
%% version to run almost everyhting
%% To plot or not to plot?
plotSCWvT=              true;  % Single Cell Width vs Time
    plotSC_fitRes=      false; % plot fit residuals?
plotSCLvT=              true;  % Single Cell Length vs Time
    plotAllLvT=         true;  % All cells Length vs Time
        plotAllnormLvT= false; % All cells normalized Length vs Time
plotdWvT=               false; % "constriction": dWidth vs time
plotBIChist=            false; % histogram of BIC of fits
plotSCWvTZ=             true;  % Single Cell Divisome and inner membrane diameter vs Time
plotCurv=               false; % Curvature fit (debug) (slow!)
plotCurvAll=            false; % All cells Curvature vs normalized time
plotCentPos=            false; % All cells normalized position of constriction site

% To fit or not to fit?
fitSC_Feingold=         true;  % Feingold model to single cell data
fitSC_Xiaofree=         true;  % Xiao model to single cell data
fitSC_2D=               false; % 2D model to single cell data
fitSCZ_Feingold=        true;  % Feingold model to single cell divisome diameter data
fitSCZ_Xiaofree=        true;  % Xiao model to single cell divisome diameter data
fitSCmts_Feingold=      true;  % Feingold model to single cell diameter data
fitSCmts_Xiaofree=      true;  % Xiao model to single cell diameter data
fitSC_Xiaofree_Tz=      false; % Xiao model with fixed Tc to single cell data
fitSC_Xiaofree_Tc_Tg =  true;  % Xiao model with fixed Tc and Tg to single cell data % Overwrites fitSC_Xiaofree_Tz
fitSC_Feingold_Tz=      false; % Feingold model with fixed Tc to single cell data
fitSC_XiaoAlpha =       false; % Xiao model with fixed Tc to single cell data

% To calculate or not to calculate?
calcCurv=               true; % Radius of curvature of cell center line
calcMTSmid=             true; % Position of central diameter local minimum (on long axis)

% ensure dependencies
if ~calcMTSmid
    plotCurv=           false;
    plotCurvAll=        false;
end
if ~calcCurv
    plotCurv=           false;
    plotCurvAll=        false;
end