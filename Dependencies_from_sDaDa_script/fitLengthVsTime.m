%-------------------------------------------------------------------------%
%                main function: fitLengthVsTime                           %
%-------------------------------------------------------------------------%
% This function allows to compute the fit of the bacteria length  versus
% time 
%
% The input of the function are:
% - l (length of the cell)
% - t (time)
% The output of the function are the fitting parameters
% - kL: decay constant 
% - l0: length at birth


function [kL, l0, fitOut, resnorm, residuals]  = fitLengthVsTime(t, l)

kLg = 0.005; % initial guess for kL [1/min]
l0g = l(1); % initial guess for length at birth [nm]

%Ag = 2000;
% p0 = [kLg, l0g, Ag];
p0 = [kLg, l0g];

% Lower bound for the parameters p0
pLB =[0, 0.5*l(1)];

% Upper bounds
pUB =[2, 1.3*l(1)];

% Solve nonlinear curve-fitting (data-fitting) problems in least-squares sense
[p, resnorm, residuals] = lsqcurvefit(@growthModel, p0, t, l, pLB, pUB);

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


%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=linear constriction
%3. (t>tg) w = 0;

% wOut = 0.*t;
% 
% %1. (t<tc), w=wMax
% isPreDiv = t<0;
% wOut(isPreDiv) = wMax;
% %wOut(isPreDiv)=1;
% 
% %2. (tc<t<tg),w=linear constriction
% %passes through points
% %(tc,wMax), (tg,0);
% isDiv = (t>=tc & t<tg);
% tDiv = t(isDiv);
% wOut(isDiv) = wMax.*sqrt( 1 - 1/(tg-tc).*(tDiv-tc) );
% 
% 
% 
% %3. (t>tg) w = 0;
% isPostDiv= (t>=tg);
% wOut(isPostDiv)=0;
end

