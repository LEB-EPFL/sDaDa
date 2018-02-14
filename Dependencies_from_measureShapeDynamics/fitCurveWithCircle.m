% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Anna Archetti
% 
% Contact:
% e-mail:
% anna.archetti@epfl.ch
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
%                main function: fitCurveWithCircleBound                         %
%-------------------------------------------------------------------------%
% This function allows to compute the curvature of a curve by fitting it
% with a circle
% 
%
% The input of the function are:
% - x, y  coordinates of the curve
%
% The output of the function are the fitting parameters
% - the radiuos of the osculating circle
% - the curvature k

function [R, xC, yC, sigmaR]  = fitCurveWithCircle(x, y, xCguess, yCguess, Rguess)

    R0 = Rguess;
    xC0 = xCguess;
    yC0 = yCguess;
    a0 = [R0, xC0, yC0];
    lb = [0, 0, 0];
    ub = [8000, Inf, Inf];

    func = @(a) circleCurve(a, x, y);
    
    [a, resnorm, residual, exitflag, output, lambda, J] = lsqnonlin(func, a0, lb, ub);
   % a = lsqcurvefit(func, a0, lb, ub);
  
   % Confidence interval 95% = 2Sigma
    aInt = nlparci(a, residual, 'jacobian', J);
    sigmaR = (aInt(1, 2) - a(1, 1))/4;
    
    R = a(1);
    xC = a(2);
    yC = a(3);

end

%----------------------
function fOut = circleCurve(a, x, y)

    R = a(1);
    xC = a(2);
    yC = a(3);

    fOut = zeros( 1, length(x));
    for i = 1:length(x)
        fOut(i) = abs((xC - x(i))^2 + (yC - y(i))^2 - R^2);
    end


end

