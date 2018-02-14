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

%-------------------------------------------------------------------------%
%                main function: fitNormLengthVsTime                           %
%-------------------------------------------------------------------------%
% This function allows to compute the fit of the bacteria length  versus
% time 
%
% The input of the function are:
% - l (length of the cell)
% - t (time)
% The output of the function are the fitting parameters
% - kL: dacay constant 
% - l0: length at birth


function [kL, l0, fitOut]  = fitNormLengthVsTime(t, wl)

kLg = 0.008; % initial guess for kL [1/min]
wl0g = wl(1); % initial guess for length at birth [nm]

%Ag = 2000;
% p0 = [kLg, l0g, Ag];
p0 = [kLg, wl0g];

% Lower bound for the parameters p0
% pLB =[0, 1, 1000];
pLB =[0, 0];

% Upper bounds
% pUB =[20, 40, 3000];
pUB =[2, 1];

% Options definition
options = optimset('Display','iter-detailed');

% Solve nonlinear curve-fitting (data-fitting) problems in least-squares sense
p = lsqcurvefit(@growthModel, p0, t, wl, pLB, pUB, options);

kL = p(1);
l0 = p(2);

tFit = unique(sort(t));
lFit = growthModel(p, tFit);
fitOut = [tFit(:), lFit(:)];
end

%----------------------
function lOut = growthModel(p, t)

kL = p(1);
l0 = p(2);

%lOut = l0.*exp(-kL.*t) + A;
lOut = l0.*exp(kL.*t);

end

