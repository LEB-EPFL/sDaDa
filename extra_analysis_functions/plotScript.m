% Copyright (C) 2017 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
% Laboratory of Experimental Biophysics
% 
% Authors: Anna Archetti and Aster Vanhecke 
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

% % This is a script for plotting results from the SIM analysis
% First manually put Tvar, Lvar, twMatrix and tlMatrix from different
% experiments together, then execute this script.
% Tvar=[Tvar30 ; Tvar31; Tvar32; Tvar33; Tvar36; Tvar37; Tvar40; Tvar42];
% Lvar=[Lvar30 ; Lvar31 ; Lvar32 ; Lvar33 ; Lvar36 ; Lvar37 ; Lvar40 ; Lvar42];
% tlMatrixAll=[tlMatrix30, tlMatrix31, tlMatrix32, tlMatrix33, tlMatrix36, tlMatrix37, tlMatrix40, tlMatrix42];
% twMatrixAll=[twMatrix30, twMatrix31, twMatrix32, twMatrix33, twMatrix36, twMatrix37, twMatrix40, twMatrix42];
% Tvar_=Tvar;
% Tvar_(11,:)=[NaN NaN NaN NaN];
% Tvar_(12,2)=NaN;

%%
% code to put all cellInfo and CellDivVar together
for frIdx=75:-1:1 %go backwards for preallocation
   cellInfo{frIdx,1}=[cellInfo07{frIdx,1} cellInfo08{frIdx,1} cellInfo10{frIdx,1} cellInfo11{frIdx,1} cellInfo13{frIdx,1} cellInfo15{frIdx,1}];
end

% change cell numbers 
divededVarEachFrm08(:,2)=divededVarEachFrm08(:,2)+max(divededVarEachFrm07(:,2));
divededVarEachFrm10(:,2)=divededVarEachFrm10(:,2)+max(divededVarEachFrm08(:,2));
divededVarEachFrm11(:,2)=divededVarEachFrm11(:,2)+max(divededVarEachFrm10(:,2));
divededVarEachFrm13(:,2)=divededVarEachFrm13(:,2)+max(divededVarEachFrm11(:,2));
divededVarEachFrm15(:,2)=divededVarEachFrm15(:,2)+max(divededVarEachFrm13(:,2));

divededVarEachFrm=[divededVarEachFrm07 ;divededVarEachFrm08 ; divededVarEachFrm10 ;divededVarEachFrm11 ;divededVarEachFrm13 ;divededVarEachFrm15];

%% Making the variable names easier:

LcF=Lvar(:,1);
LgF=Lvar(:,2);
LcX=Lvar(:,4);
LgX=Lvar(:,5);
Lc2D=Lvar(:,6);
Lg2D=Lvar(:,7);


TcF=Tvar(:,1);
TgF=Tvar(:,2);
TcX=Tvar(:,3);
TgX=Tvar(:,4);
Tc2D=Tvar(:,5);
Tg2D=Tvar(:,6);

Lb=Lvar(:,3);

try
    Lz=Lvar(:,8);
    Tz=Tvar(:,7);
catch
end

% Let's use the Feingold variables as the defaults
Lc=LcF;
Lg=LgF;
Tc=TcF;
Tg=TgF;

%%
% max(divededVarEachFrm30(:,2))+max(divededVarEachFrm31(:,2))+max(divededVarEachFrm32(:,2))+max(divededVarEachFrm33(:,2))+max(divededVarEachFrm36(:,2))+max(divededVarEachFrm37(:,2))+max(divededVarEachFrm40(:,2))+max(divededVarEachFrm42(:,2));

nbins=10;
% close all

%Histograms
figure('Name','Times Histogram','NumberTitle','off'),
subplot(2,3,1), hist(Tvar,nbins), legend('Tc 2D','Tg 2D','Tc Feingold','Tg Feingold')
% subplot(2,3,2), hist(Tc2D,nbins), legend('Tc 2D model Histogram')
% subplot(2,3,3), hist(Tg2D,nbins), legend('Tg 2D model Histogram')
subplot(2,3,5), hist(TcF,nbins), legend('Tc Feingold model Histogram')
subplot(2,3,6), hist(TgF,nbins), legend('Tg Feingold model Histogram')
subplot(2,3,4), hist(TgF-TcF, nbins), legend('tau Feingold')
subplot(2,3,2), hist(TcX,nbins), legend('Tc Xiao model Histogram')
subplot(2,3,3), hist(TgX,nbins), legend('Tg Xiao model Histogram')


figure('Name','Lengths Histogram','NumberTitle','off'),
subplot(2,3,1), hist(Lvar,nbins), legend('Lc 2D','Lg 2D','Lc Feingold','Lg Feingold')
% subplot(2,3,2), hist(Lc2D,nbins), legend('Lc 2D model Histogram')
% subplot(2,3,3), hist(Lg2D,nbins), legend('Lg 2D model Histogram')
subplot(2,3,5), hist(LcF,nbins), legend('Lc Feingold model Histogram')
subplot(2,3,6), hist(LgF,nbins), legend('Lg Feingold model Histogram')
subplot(2,3,2), hist(LcX,nbins), legend('Lc Xiao model Histogram')
subplot(2,3,3), hist(LgX,nbins), legend('Lg Xiao model Histogram')

%Correlations?
corrFig(Tc,Tg);
corrFig(Tg-Tc,Lc,'Lc vs Tg-Tc');
corrFig(Tg-Tc,Lg-Lc,'Lg-Lc vs Tg-Tc');
corrFig(Tg-Tc,Lg,'Lg vs Tau');
corrFig(Lc,Lg);
corrFig(Lg-Lc,Lg,'Lg vs Lg-Lc');
corrFig(Lc, Lg-Lc,'Lg-Lc vs Lc');
corrFig(Tc,Tg-Tc,'Tg-Tc vs Tc');
corrFig(Tg,Tg-Tc,'Tg-Tc vs Tg');


%%
figure('Name','Lb correlations','NumberTitle','off'),
% subplot(2,3,1), plot(Lb,Lc,'o'), legend('Lc 2D vs Lb')
subplot(2,3,2), plot(Lb,LcF,'o'), legend('Lc Feing vs Lb')
% subplot(2,3,3), plot(Lb,Lg,'o'), legend('Lg 2D vs Lb')
subplot(2,3,4), plot(Lb,LgF,'o'), legend('Lg Feing vs Lb')
% subplot(2,3,5), plot(Lb,Lg-Lc,'o'), legend('Elongation during constriction 2D vs Lb')
subplot(2,3,6), plot(Lb,LgF-LcF,'o'), legend('Elongation during constriction Feing vs Lb')

subplot(2,3,1), plot(Lb,LcX,'o'), legend('Lc Xiao vs Lb')
subplot(2,3,3), plot(Lb,LgX,'o'), legend('Lg Xiao vs Lb')
subplot(2,3,5), plot(Lb,LgX-LcX,'o'), legend('Elongation during constriction Xiao vs Lb')

figure('Name','Lb correlations with time','NumberTitle','off'),
subplot(2,3,1), plot(Lb,Tc,'o'), legend('Tc 2D vs Lb')
subplot(2,3,2), plot(Lb,TcF,'o'), legend('Tc Feing vs Lb')
subplot(2,3,3), plot(Lb,Tg,'o'), legend('Tg 2D vs Lb')
subplot(2,3,4), plot(Lb,TgF,'o'), legend('Tg Feing vs Lb')
subplot(2,3,5), plot(Lb,Tg-Tc,'o'), legend('tau vs Lb')
% subplot(2,3,6), plot(Lb,TgF-TcF,'o'), legend('tau Feing vs Lb')

corrFig(Lc-Lb, Lg-Lc,'Lg-Lc vs Lc-Lb')

% cte size extension? Lg-Lb is "size extension"

corrFig(Lb, Lg-Lb,'Lg-Lb vs Lb')

% corrFig(Lb, LgF-Lb,'Lg-Lb vs Lb Feingold')

% Size extension seems to be mostly uncorrelated with Lb, except for the
% edge cases: really small cells (1,7-1,8µm) elongate more (0.9-1.4µm), really long cells (2.2-2.4µm)do not
% elongate a lot(0.7-0.9µm).

%intergenerational continuity /inheritance of cell size
figure('Name','Lg/2 vs Lb Feingold','NumberTitle','off'),
plot(Lb, LgF/2,'o')

% Xpos=Lb;
% 
% [f3,x3] = ksdensity(Xpos, 'Bandwidth',50);
% [f4,x4] = ksdensity(Xpos, 'Bandwidth',60);
% [f5,x5] = ksdensity(Xpos, 'Bandwidth',70);
% [f6,x6] = ksdensity(Xpos, 'Bandwidth',80);
% [f7,x7] = ksdensity(Xpos, 'Bandwidth',90);
% [f8,x8] = ksdensity(Xpos, 'Bandwidth',100);


figure('Name','Lb Histogram','NumberTitle','off'), hist(Lb,nbins),
% hold on
% plot(x3,6000*f3,'r');
% plot(x4,6000*f4,'g');
% plot(x5,6000*f5,'k');
% plot(x6,6000*f6,'c');
% plot(x7,6000*f7,'m');
% plot(x8,6000*f8,'b');
% legend('hist','50', '60', '70', '80', '90', '100')
% hold off

%% Execute other plot scripts
PlotAllCellsWvsT
PlotFit10cells

%% Plotting Elongation rate and constriction rate
figure, plot(log(2)./(ElMat(:,1)),Tg-Tc,'o')
title('Duration of constriction vs elongation rate', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Elongation doubling time (min)','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Duration of constriction (min)','FontSize',16, 'FontName', 'Calibri Light' );
% xlim([100 300]);

 %% Jie Xiao's plots
% 
% % Vc/Vec
% 
% figure;
% plot((Lg-Lc)./(Tg-Tc),500./(Tg-Tc),'o')
% title('Vc/Vec', 'FontSize', 20, 'FontName', 'Calibri Light'),
% xlabel('Elongation rate after constriction (nm/min)','FontSize',16, 'FontName', 'Calibri Light'  );
% ylabel('Constriction rate (nm/min)','FontSize',16, 'FontName', 'Calibri Light' );
% % xlim([0 15]);
% % ylim([0 15]);
% 
% VcoverVec=(500./(Tg-Tc))./((Lg-Lc)./(Tg-Tc));
% figure,
% boxplot(VcoverVec)
% % Vc/Vep
% figure;
% plot((Lc-Lb)./Tc,500./(Tc),'o')
% title('Vc/Vep', 'FontSize', 20, 'FontName', 'Calibri Light'),
% xlabel('Elongation rate before constriction (nm/min)','FontSize',16, 'FontName', 'Calibri Light'  );
% ylabel('Constriction rate (nm/min)','FontSize',16, 'FontName', 'Calibri Light' );
% % xlim([0 35]);
% % ylim([0 35]);
% hold all;
% 
% VcoverVep=(500./(Tg-Tc))./((Lc-Lb)./(Tc));
% figure,
% boxplot(VcoverVep)
 
 %% Vep vs Vec
% figure,
% plot(((Lc-ElMat(:,2))./(Tc)),(Lg-Lc)./(Tg-Tc),'o')
% xlim([0 25]);
% ylim([0 15]);
% title('Vep vs Vec', 'FontSize', 20, 'FontName', 'Calibri Light'),
