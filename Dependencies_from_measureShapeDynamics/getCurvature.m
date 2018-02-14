function [k, dfdx, d2fdx2, x, y] = getCurvature(xTemp, yTemp)

%     % Grab the angle that the tg at the central point form with the x axis  
%     dX = xTempOrig( round(length(xTempOrig)/2)) - xTempOrig( round(length(xTempOrig)/2)+1 ) ;
%     dY = yTempOrig( round(length(yTempOrig)/2)) - yTempOrig( round(length(yTempOrig)/2)+1 );
%     
%     theta = atan2(dY,dX);
% 
%     % Rotate parallel to xAxis 
%     thetaRot = -theta + pi/3;
%     R = [cos(thetaRot),-sin(thetaRot );...
%     sin(thetaRot ),cos(thetaRot)];
% 
%     xTemp = xTempOrig.*cos(thetaRot) - yTempOrig.*sin(thetaRot);
%     yTemp = xTempOrig.*sin(thetaRot) + yTempOrig.*cos(thetaRot);
    
    % Debug
%     figure, 
%     plot(xTempOrig, yTempOrig, 'm*'), hold on
%     plot(xTemp, yTemp, 'k*'),  axis equal
%     legend('Original', 'rotated')
    
    % Make the interpolation of the central part of the bacteria central
    % line
    xTempTemp = xTemp(length(xTemp)/2 - length(xTemp)/4:length(xTemp)/2 + length(xTemp)/4);
    yTempTemp = yTemp(length(yTemp)/2 - length(yTemp)/4:length(yTemp)/2 + length(yTemp)/4);

    xyTemp = sortrows([xTempTemp yTempTemp], 1);
    xyTemp = sortrows(xyTemp, 2);

    xTempForCheck = xyTemp(:,1);
    yTempForCheck = xyTemp(:,2);
    % obs: xTemp is ordered in a crescent order respect to x and then respect
    % to y

    % DEBUG
    % check if monothonic
%     figure, plot(xTemp, yTemp, 'm*'),  axis equal
%     figure, plot(xTempTemp, yTempTemp, 'go'),  axis equal
    
    % Check if there are x_n+1 value smaller than x_n and the corresponding
    % y_n+1 value bigger than y_n
%     xDiff = diff(xTempForCheck);
%     yDiff = diff(yTempForCheck);
    xDiff = diff(xTempTemp);
    yDiff = diff(yTempTemp);
    
    % if the differences in y are poistive and all the differences in x are
    % negative it means that the function is monotonic decrescent and the
    % data points are flipped
    if size(find(xDiff <= 0)) == size(xDiff) 
            
        xTempTemp = fliplr(xTempTemp);
        yTempTemp = fliplr(yTempTemp);
    end
    monCheckX = find(xDiff <= 0);
    monCheckY = find(yDiff(monCheckX) >= 0);
    %if nnz(monCheckX) >0 && nnz(monCheckY) >0 && nnz(monCheckX)<length(xTempForCheck)-2
     if nnz(monCheckX) >0 
       %figure, plot(xTemp, yTemp, 'm*'),  axis equal
       xForRotOrig = xTemp;
       yForRotOrig = yTemp;
       xForRot = xTempTemp;
       yForRot = yTempTemp;
       
       xTemp = yForRotOrig;
       yTemp = -xForRotOrig;
       xTempTemp  = yForRot;
       yTempTemp  = -xForRot;

       % DEBUG
%        figure, plot(xTemp, yTemp, 'r*'),  axis equal
%        figure, plot(xTempTemp, yTempTemp, 'go'),  axis equal
     end

     if isempty( find(diff(yTempTemp)<0 | diff(xTempTemp)<0 ,1)) % if the function is monothonic
         % TODO: check if it is monotonically INCREASING
         fitMod = griddedInterpolant(xTempTemp, yTempTemp, 'spline');
         % Not more 5nm precision expected
         %x = (min(xTempTemp):20:max(xTempTemp));
         xx = xTempTemp;
         yy = fitMod(xx);
         
         y=smooth(xx, yy, 'rloess');
         [x,ind] = sort(xx);
         y = y(ind);
         % Compute the first and second derivative
         [dfdx, d2fdx2] = dfdxc(x,y');
         %k = sqrt(d2fdx2.^2);
         
         % Compute the curvature k
         k = sqrt(d2fdx2.^2).*(1./(1 + (dfdx).^2).^(3/2));
     else
         
         xx = xTempTemp;
         yy = yTempTemp;
         y=smooth(xx, yy, 'rloess');
         [x,ind] = sort(xx);
         y = y(ind);
         % Compute the first and second derivative
         [dfdx, d2fdx2] = dfdxc(x,y');
         %k = sqrt(d2fdx2.^2);
         
         % Compute the curvature k
         k = sqrt(d2fdx2.^2).*(1./(1 + (dfdx).^2).^(3/2));
     end
    
    
    % Set figure properties
%     scrsz = get(0,'ScreenSize');
%     h = figure('Position',[(scrsz(3)-1280)/3 (scrsz(4)-720)/3 1280 720],'color','w');
%     xCirc = (0:pi/64:2*pi);
    % Plot the osculating circle tangent to the curve
    %for pointTg = 2: length(x)-1
%-------------------------------------------------------------------------%
% DEBUG
%     % Plot the osculationg circle corresponding to the tagent point with
%     % highest curvature
%     pointTg =  round(length(k)/2) - round(length(k)/4) + find(max( k(round(length(k)/2) - round(length(k)/4):round(length(k)/2) + round(length(k)/4)) ));
%     set(0,'CurrentFigure',h);
% 
%     % Smoothed line
%     plot(x, y, 'b--', 'LineWidth', 5), 
%     axis equal
%     hold on
%     
%     % Tangent point
%     plot(x(pointTg), y(pointTg), 'k*', 'MarkerSize',10)
%     
%     % Interpolation line
%     plot(xx, yy, 'g.-'),
%     
%     % Septum points
%     plot(xTempTemp, yTempTemp, 'ro')
%     
%     % All edge points
%     plot(xTemp, yTemp, 'mo')
%     % Plot osculating circle
%     R = 1/k(pointTg);
%     if d2fdx2(pointTg) > 0
%        if  dfdx(pointTg) > 0
%             Xc = x(pointTg) - R*sin( atan(dfdx(pointTg)) );
%             Yc = y(pointTg) + R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%        else
%             Xc = x(pointTg) - R*sin( atan(dfdx(pointTg)) );
%             Yc = y(pointTg) + R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%        end
%         
%     else
%         if dfdx(pointTg) > 0
%             Xc = x(pointTg) + R*sin( atan(dfdx(pointTg)) );
%             Yc = y(pointTg) - R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%         else
%             Xc = x(pointTg) + R*sin( atan(dfdx(pointTg)) );
%             Yc = y(pointTg) - R*cos( atan(dfdx(pointTg)) );
%             plot(R*cos(xCirc) + Xc, R*sin(xCirc) + Yc, 'r--')
%             axis equal
%             
%         end
%     end
    % END DEBUG
%-------------------------------------------------------------------------%
%     pointTg
%     k(pointTg)
%    drawnow;
% end
end