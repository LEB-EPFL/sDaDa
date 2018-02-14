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

function [cumSurface,cumVolume] = compute3DShape2(X,Y)
% More general take on calculating surface area and volume of a radially
% symetric shape described by vextors X and Y, with Y being the radius at
% different lengths in vector X. X and Y should be vectors of equal length.
% The shape is subdivided in conical frusta, to calculate surface area's
% and volumes of the subdivisions, of which the cumulative sum is
% calculated.

% test input
narginchk(2,2)
if size(X) ~= size(Y)
    error('dimensions of X and Y do not match!')
end

X=X(~isnan(Y));
Y=Y(~isnan(Y));

% For code readability:
h  = (X(2:end)-X(1:end-1));
R1 = Y(1:end-1);
R2 = Y(2:end);

Surface=pi.*(R1+R2).*sqrt( (R1-R2).^2 + h.^2); % Surface area of conical frustum
Volume=pi/3.*h.*(R1.^2 + R1.*R2 + R2.^2); % Volume of conical frustum

cumSurface=cumsum(Surface);
cumVolume=cumsum(Volume);
end