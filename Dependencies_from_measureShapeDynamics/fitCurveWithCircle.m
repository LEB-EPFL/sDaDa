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

