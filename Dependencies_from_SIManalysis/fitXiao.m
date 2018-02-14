function [tc, tg, wMax, alpha, fitOut, resnorm, residual, BIC]  = fitXiao(t,w,varargin)
% function [tc, tg, wMax, alpha, fitOut, resnorm, residual, BIC]  = fitXiao(t, w)
% This function allows to compute the fit of the waist diameter, w, versus
% time, t, to the Jiao model Xiao
% The input of the function are:
% - t: vector timepoints
% - w: vector containing diameter at septum divided by the max diameter for
% the respective timepoints
% Optional inputs:
% fitXiao(t,w,'alpha',a): fixes the alpha parameter to the value of a
% fitXiao(t,w,'tc',tc): fixes the tc parameter to the value to tc
% The output of the function are the fitting parameters, tc, tg, wMax and
% alpha, as weel as the points making up the fitted curve. Optional output:
% resnorm and residual of the lsqcurvefit and the Bayesian information
% criterion.
narginchk(2,6);

nparam=4;
lAlpha=0.1;
uAlpha=10;
alpha0 = 2;
tc0 = 0;
ltc=0;
utc=max(t);
tg0 = max(t);
ltg= 0;
utg = Inf;
if size(varargin,2) ~= 0
    varIn=1;
    while varIn<numel(varargin)
        switch varargin{1,varIn}
            case {'alpha', 'Alpha'}
                lAlpha=varargin{1,varIn+1};
                uAlpha=lAlpha;
                nparam=nparam-1;
                alpha0 = lAlpha;
                varIn=varIn+2;
            case {'tc','Tc'}
                nparam=nparam-1;
                ltc=varargin{1,varIn+1};
                utc=ltc;
                tc0=ltc;
                varIn=varIn+2;
            case {'tg','Tg'}
                nparam=nparam-1;
                ltg=varargin{1,varIn+1};
                utg=ltg;
                tg0=ltg;
                varIn=varIn+2;
            otherwise
                error('Invalid input')
        end
    end
end
if w(1)<2
    wMax0=1;
    lwMax0=0.7;
    uwMax0=1.5;
else
    wMax0=500;
    lwMax0=350;
    uwMax0=800;
end

a0 = [tc0,tg0,wMax0, alpha0];
lb =[ltc,ltg,lwMax0,lAlpha];
ub =[utc,utg,uwMax0,uAlpha];
[a, resnorm, residual] = lsqcurvefit(@XiaoModel, a0, t, w, lb, ub);
tc = a(1);
tg = a(2);
wMax = a(3);
alpha = a(4);

% Compute the Bayesian information criterion
numDataPoints = length(residual);
BIC = numDataPoints*log(resnorm)+nparam*log(numDataPoints);

tFit = unique(sort(t));
wFit = XiaoModel(a,tFit);
fitOut = [tFit(:),wFit(:)];
end

%----------------------
function wOut = XiaoModel(a, t)

tc = a(1);
tg = a(2);
wMax = a(3);
alpha = a(4);
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
wOut(isDiv) = wMax.*( 1 - ((tDiv-tc)./(tg-tc)).^alpha ).^(1./alpha);


%3. (t>tg) w = 0;
isPostDiv= (t>=tg);
wOut(isPostDiv)=0;

%4. tc<tg, so punish when tg<t<tc
isYouDoneFdUp = (t>=tg & t<tc);
wOut(isYouDoneFdUp)=-9001;
end