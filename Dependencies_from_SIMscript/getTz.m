function Tz=getTz(ZassemVec,tVec) % ZassemVec=tZmat(:,cellIdx*4), tVec=tZmat(:,(cellIdx*4)-3)
if sum(ZassemVec)
    SE=strel('line',3,90);
    ZassemVec=imclose(ZassemVec,SE);
    ZassemVec=imopen(ZassemVec,SE); % remove single positives
    % imclose(ZassemVec,SE); %fill holes
    frIdx=find(ZassemVec,1);
    Tz=tVec(frIdx);
else
    Tz=NaN;
end
if isempty(Tz)
    Tz=NaN;
end
end
