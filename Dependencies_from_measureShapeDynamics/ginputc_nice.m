% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Seamus Holden, Anna Archetti and Aster Vanhecke 
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

function [xOut,yOut,button] = ginputc_nice(imsize,varargin)
%get rid of points outside of image

[x,y,button] = ginputc(varargin{:});

nPt = numel(x);
xOut =[];
yOut =[];
buttonOut =[];
for ii=1:nPt
  if y(ii)>=1 && y(ii)<=imsize(1) && x(ii)>=1 && x(ii)<=imsize(2)
    xOut=[xOut,x(ii)];
    yOut=[yOut,y(ii)];
    buttonOut=[buttonOut,button(ii)];
  end
end
