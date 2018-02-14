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

function Tz=getTz(ZassemVec,tVec) % ZassemVec=tZmat(:,cellIdx*4), tVec=tZmat(:,(cellIdx*4)-3)
if sum(ZassemVec)
    SE=strel('line',3,90);
    ZassemVec=imclose(ZassemVec,SE);
    ZassemVec=imopen(ZassemVec,SE); % remove single positives
    % imclose(ZassemVec,SE); %fill holes
    frIdx=find(ZassemVec,1);
    Tz=tVec(frIdx);
else
    Tz=NaN;
end
if isempty(Tz)
    Tz=NaN;
end
end
