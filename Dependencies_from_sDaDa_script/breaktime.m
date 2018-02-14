% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Aster Vanhecke (based on a MATLAB tutorial)
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
function [weeks, days, hours, mins, secs]=breaktime(totseconds, varargin)
%breaktime breaks a number of seconds into hours, minutes and seconds
%Format: [weeks, days, hours, mins, secs]= breaktime(totalseconds)

weeks=floor(totseconds/604800);
remsecs=rem(totseconds,604800);
days=floor(remsecs/86400);
remsecs=rem(remsecs,86400);
hours=floor(remsecs/3600);
remsecs=rem(remsecs,3600);
mins=floor(remsecs/60);
secs=rem(remsecs,60);
if ~isempty('varargin')
    if strcmp(varargin, 'disp')
        disp([num2str(weeks) ' weeks, ' num2str(days) ' days, ' num2str(hours) ' hours, ' num2str(mins) ' minutes and ' num2str(secs) ' seconds.'])
    end
end
end