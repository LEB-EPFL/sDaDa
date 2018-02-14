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


