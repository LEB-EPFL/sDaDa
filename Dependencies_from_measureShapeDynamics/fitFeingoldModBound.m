%-------------------------------------------------------------------------%
%                main function: fitFeingoldModBound                         %
%-------------------------------------------------------------------------%
% This function allows to compute the fit of the waist diamtert w versus
% time considering a model where
% - the area inserted at septum is a 3D hemisferic
% - the area rate insertion is constant
%
% The input of the function are:
% - w (diameter at septum divided by the max diameter)
% - t (time)
% - optional: tc
% The output of the function are the fitting parameters

function [tc, tg, wMax, fitOut, resnorm, residual, BIC]  = fitFeingoldModBound(t,w,varargin)
narginchk(2,6);
nparam=3;
if size(varargin,2)== 0
    tc0 = min(t);
    ltc=0;
    utc=max(t);
else
    tc0=varargin{1,1};
    ltc=tc0;
    utc=tc0;
    nparam=nparam-1;
end
if w(1)<2
    wMax0=1;
    lwMax0=0.9;
    uwMax0=1.2;
else
    wMax0=500;
    lwMax0=350;
    uwMax0=800;
end
tg0 = max(t);
a0 = [tc0,tg0,wMax0];
lb =[ltc,0,lwMax0];
ub =[utc,Inf,uwMax0];
[a, resnorm, residual] = lsqcurvefit(@FeingoldModel, a0, t, w, lb, ub);
tc = a(1);
tg = a(2);
wMax = a(3);

% Compute the Bayesian information criterion
numDataPoints = length(residual);
BIC = numDataPoints*log(resnorm)+nparam*log(numDataPoints);

tFit = unique(sort(t));
wFit = FeingoldModel(a,tFit);
fitOut = [tFit(:),wFit(:)];
end

%----------------------
function wOut = FeingoldModel(a, t)

tc = a(1);
tg = a(2);
wMax = a(3);
%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=linear constriction
%3. (t>tg) w = 0;
%4. tc<tg, so punish when tg<t<tc

wOut = 0.*t;

%1. (t<tc), w=wMax
isPreDiv = t<tc;
wOut(isPreDiv) = wMax;
%wOut(isPreDiv)=1;

%2. (tc<t<tg),w=linear constriction
%passes through points
%(tc,wMax), (tg,0);
isDiv = (t>=tc & t<tg);
tDiv = t(isDiv);
wOut(isDiv) = wMax.*sqrt( 1 - (tDiv-tc).^2 ./(tg-tc).^2 );



%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;

%4. tc<tg, so punish when tg<t<tc
isYouDoneFdUp = (t>=tg & t<tc);
wOut(isYouDoneFdUp)=-9001;
end

