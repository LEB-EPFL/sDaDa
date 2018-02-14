%-------------------------------------------------------------------------%
%                main function: fitSeamusModBound                         %
%-------------------------------------------------------------------------%
% This function allows to compute the fit of the waist diamtert w versus
% time considering a model where
% - the area inserted at septum is a 2D circle
% - the area rate insertion is constant
%
% The input of the function are:
% - w (diameter at septum divided by the max diameter)
% - t (time)
% The output of the function are the fitting parameters

function [tc, tg, wMax, fitOut, resnorm, residual]  = fitSeamusModBound(t, w)

tc0 = 0;
tg0 = max(t);
wMax0 = 1;
a0 = [tc0,tg0,wMax0];
lb =[0,0,0];
ub =[max(t),Inf,1];
[a, resnorm, residual]= lsqcurvefit(@seamusModel, a0, t, w, lb, ub);
tc = a(1);
tg = a(2);
wMax = a(3);

tFit = unique(sort(t));
wFit = seamusModel(a,tFit);
fitOut = [tFit(:),wFit(:)];
end

%----------------------
function wOut = seamusModel(a, t)

tc = a(1);
tg = a(2);
wMax = a(3);
%fic

%3 conditions
%1. (t<tc), w=wMax
%2. (tc<t<tg), w=linear constriction
%3. (t>tg) w = 0;

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
wOut(isDiv) = wMax.*sqrt( 1 - 1/(tg-tc).*(tDiv-tc) );



%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;
end

