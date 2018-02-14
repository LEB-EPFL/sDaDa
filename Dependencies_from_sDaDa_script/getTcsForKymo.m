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

% Script to prepare a matrix with Tc's for plotting on the kymographs. It
% takes the Tc values for Feingold and Xiao either from Tvar or from
% fitResultX/F, also gets Tz from fitResultFz or fitResultXz
TcsForKymo=[];
if exist('Tvar','var')
    TcsForKymo=Tvar(:,1);
    TcsForKymo(:,2)=Tvar(:,3);
else
    if exist('fitResultF','var')
        TcsForKymo(:,end+1)=fitResultF(1,:)';
    end
    if exist('fitResultX','var')
        TcsForKymo(:,end+1)=fitResultX(1,:)';
    end
end
try
    TcsForKymo(:,end+1)=fitResultFz(1,:)';
catch
    TcsForKymo(:,end+1)=fitResultXz(1,:)';
end