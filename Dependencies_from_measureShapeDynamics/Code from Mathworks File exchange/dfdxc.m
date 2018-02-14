% -------------------------------------------------------------- dfdxc(x,f)
% DFDXC(x,f) computes the first two derivatives of f(x) using a modified
%            form of the central difference formula.
%
% Inputs:
%
%       x = [I,1] or [1,I] vector. independent variable.
%       f = [I,1] or [1,I] vector.   dependent variable, f(x)
%
% Outputs:
%
%       dfdx   ==   d(f(x))/dx      first  derivative
%       d2fdx2 == d^2(f(x))/dx^2    second derivative
%
% Example:
%       x = [0 : 0.1 : 1]';
%       f = x.^3;
%       
%       [dfdx, d2fdx2] = dfdxc(x,f);
% 
% 
%       figure; hold on, grid on,  plot(x,dfdx  ,'r'),  plot(x,3*x.^2,'k.'), 
%       figure; hold on, grid on,  plot(x,d2fdx2,'r'),  plot(x,6*x   ,'k.'), 
%
%
% See also:
%   DIFF, PCHIP.
% -------------------------------------------------------------------------

function [dfdx, d2fdx2] = dfdxc(x,f)

[Ix,Jx] = size(x);
[If,Jf] = size(f);

if (Ix ~= If) || (Jx ~= Jf)
    disp('ERROR: f(x) must be the same size as x.') 
    return,
end

if (Ix ~= 1) & (Jx ~= 1)
    disp('ERROR: x must be a vector.')
    return, 
end


% ------------------------ Central divided difference formula for the slope
dfdx   = pchip(0.5*(x(1:end-1)+x(2:end)), diff(f)   ./diff(x), x);

% -------------------- Central divided difference formula for the curvature
d2fdx2 = pchip(0.5*(x(1:end-1)+x(2:end)), diff(dfdx)./diff(x), x);



