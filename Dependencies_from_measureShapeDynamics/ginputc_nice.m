function [xOut,yOut,button] = ginputc_nice(imsize,varargin)
%get rid of points outside of image

[x,y,button] = ginputc(varargin{:});

nPt = numel(x);
xOut =[];
yOut =[];
buttonOut =[];
for ii=1:nPt
  if y(ii)>=1 && y(ii)<=imsize(1) && x(ii)>=1 && x(ii)<=imsize(2)
    xOut=[xOut,x(ii)];
    yOut=[yOut,y(ii)];
    buttonOut=[buttonOut,button(ii)];
  end
end
