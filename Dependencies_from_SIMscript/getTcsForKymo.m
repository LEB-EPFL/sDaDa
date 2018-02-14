% Script to prepare a matrix with Tc's for plotting on the kymographs. It
% takes the Tc values for Feingold and Xiao either from Tvar or from
% fitResultX/F, also gets Tz from fitResultFz or fitResultXz
TcsForKymo=[];
if exist('Tvar','var')
    TcsForKymo=Tvar(:,1);
    TcsForKymo(:,2)=Tvar(:,3);
else
    if exist('fitResultF','var')
        TcsForKymo(:,end+1)=fitResultF(1,:)';
    end
    if exist('fitResultX','var')
        TcsForKymo(:,end+1)=fitResultX(1,:)';
    end
end
try
    TcsForKymo(:,end+1)=fitResultFz(1,:)';
catch
    TcsForKymo(:,end+1)=fitResultXz(1,:)';
end