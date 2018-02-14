function [cumSurface,cumVolume] = compute3DShape2(X,Y)
% More general take on calculating surface area and volume of a radially
% symetric shape described by vextors X and Y, with Y being the radius at
% different lengths in vector X. X and Y should be vectors of equal length.
% The shape is subdivided in conical frusta, to calculate surface area's
% and volumes of the subdivisions, of which the cumulative sum is
% calculated.
% Author: Aster Vanhecke

% test input
narginchk(2,2)
if size(X) ~= size(Y)
    error('dimensions of X and Y do not match!')
end

X=X(~isnan(Y));
Y=Y(~isnan(Y));

% For code readability:
h  = (X(2:end)-X(1:end-1));
R1 = Y(1:end-1);
R2 = Y(2:end);

Surface=pi.*(R1+R2).*sqrt( (R1-R2).^2 + h.^2); % Surface area of conical frustum
Volume=pi/3.*h.*(R1.^2 + R1.*R2 + R2.^2); % Volume of conical frustum

cumSurface=cumsum(Surface);
cumVolume=cumsum(Volume);
end