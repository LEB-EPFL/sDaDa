% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Seamus Holden, Anna Archetti and Aster Vanhecke 
% 
% Contact:
% e-mail:
% anna.archetti@epfl.ch
% aster.vanhecke@epfl.ch
% 
% paper mail:
% EPFL SB IPHYS LEB 
% BSP 428 (Cubotron UNIL) 
% Rte de la Sorge 
% CH-1015 Lausanne
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function res0 = watershedSeedSeg4(res,minCellSize)

[x,y] = ginputc_nice(size(res),'ShowPoints',true,'ConnectPoints',false,'Color','b');
xI = round(x);
yI = round(y);

labelIm = bwlabel(res);
zeroIm = logical(zeros(size(res)));
W = zeroIm;
seed = zeroIm;
nPt = numel(xI);
for ii =1:nPt
  cellNo = labelIm(yI(ii),xI(ii));
  W(labelIm == cellNo) = 1;
  seed(yI(ii),xI(ii))=1;
end

%watershed to split cells
gc=~W;
D=bwdist(gc);
D2 = imimposemin(-D,seed);
L=watershed(D2);
s=L==0;
imCellSplit = W&~s;
%delete small regions that may have appeared
%remove isolated bg points
imCellSplit=imfill(imCellSplit,'holes');
%remove regions too small to be cells
imCellSplit = bwareaopen(imCellSplit,minCellSize);

%update the all cell image
res0 = res;
res0(W==1 & imCellSplit==0)=0;


