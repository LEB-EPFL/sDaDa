function [width, xLeft, xRight,yLeft,yRight,xLeftMax,xRightMax]= getFWXM(varargin)
% Function [width, xLeft, xRight,yLeft,yRight,xLeftMax,xRightMax]=
% getFWXM(x,y, 'options')
% finds full width at X* the max (X=0.5 for full width at half max), 
% for (x,y) data than can have multiple peaks.
% options:
% - 'Half', X: Instead of FWHM, find FW at X*max
% - 'nbreaks', n : number of breakpoints for splinefitting to find
% maxima.
% - 'searchStep': step size for search for FWXM
% - 'FigZ' 


% Input: x,y,searchstep, half
narginchk(2,10);
x=varargin{1,1};
y=varargin{1,2};
nbreaks=10;
Half=0.5;
searchStep=0.1;
% loop to check and asign input variables
for i=3:2:size(varargin,2)-1
    switch varargin{1,i}
        case 'Half'
            Half=varargin{1,i+1};
        case 'nbreaks'
            nbreaks=varargin{1,i+1};
        case 'searchStep'
            searchStep=varargin{1,i+1};
        case 'FigZ'
            fig=varargin{1,i+1};
        otherwise
            error(['invalid input: ' varargin{1,i} ])
    end
end

xgrid=ndgrid(x(1):searchStep:x(end)); % create grid
pp=splinefit(x,y,nbreaks);
spline=ppval(pp,xgrid); % We could replace the splinefit by a gridded
% interpolant, but the data generally has this step behaviour which we want
% to get rid of with the spline.

[maxI,maxX]=findpeaks(spline,'MinPeakHeight',max(spline)/3); % maxX holds 
% id's for maxima, not position values!

% from found peaks, filter out the real left and right peaks
switch length(maxX)
    case 0 % what if no peaks are found? The maxima are at the edges.
        width=NaN;
        xLeft=NaN;
        xRight=NaN;
        yLeft=NaN;
        yRight=NaN;
        xLeftMax=NaN;
        xRightMax=NaN;
        return
    case 1 % only one peak
        posLeftMax=maxX;
        valLeftMax=maxI;
        posRightMax=maxX;
        valRightMax=maxI;
    case 2 % two peaks
        if maxX(1)<maxX(2) % if maxX(1) is left peak
            posLeftMax=maxX(1);
            valLeftMax=maxI(1);
            posRightMax=maxX(2);
            valRightMax=maxI(2);
        else
            posLeftMax=maxX(2);
            valLeftMax=maxI(2);
            posRightMax=maxX(1);
            valRightMax=maxI(1);
        end
    otherwise % more than two peaks...
        % sort peaks, pick most left and most right peaks.
        [maxX, idx]=sort(maxX,1,'ascend');
        maxI=maxI(idx);
        posLeftMax=maxX(1);
        valLeftMax=maxI(1);
        posRightMax=maxX(end);
        valRightMax=maxI(end);
%         % Debug
%         figure,
%         plot(x,y)
%         hold on, plot(x,spline)
%         plot(posLeftMax,valLeftMax,'+')
%         plot(posRightMax,valRightMax,'+')
%         hold off
end

% From left peak, scan left to find half max
posLeft = posLeftMax;
valCur = Inf;
while posLeft>1 && valCur > valLeftMax*Half % posLeft is element id in xgrid/spline!
    posLeft =posLeft-1;
    valCur = spline(posLeft);
end
% From right peak, scan right to find half max
posRight = posRightMax;
valCur = Inf;
while posRight<length(xgrid) && valCur > valRightMax*Half
    posRight =posRight+1;
    valCur = spline(posRight);
end

%%DEBUG
%hold off;
% plot(x,y);
% hold all;
% plot(posLeftMax,valLeftMax,'ro');
% plot(posLeft,valLeftMax/2,'bo')
% plot(posRightMax,valRightMax,'ro');
% plot(posRight,valRightMax/2,'bo');
%%DEBUG

xLeft = xgrid(posLeft);
xRight = xgrid(posRight);
width = xRight-xLeft;
yLeft =valLeftMax*Half;
yRight =valRightMax*Half;
xLeftMax=xgrid(posLeftMax);
xRightMax=xgrid(posRightMax);

% Figure for measureZ.m in the SIM analysis pipeline.
if exist('fig','var')
    figure(fig),
    plot(x,y)
    hold on,
    plot(xgrid,spline)
    plot(xLeftMax,valLeftMax,'+')
    plot(xRightMax,valRightMax,'+')
    plot(xLeft,yLeft,'+')
    plot(xRight,yRight,'+')
    ylim([0 Inf])
    hold off
    try
        title(['Width: ' num2str(width) 'pxl or ' num2str(width*30) ' nm'])
    catch
    end
end

end
