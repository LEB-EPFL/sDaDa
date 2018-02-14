function [twMatrix, tlMatrix, Tvar, Lvar, BIC, p1C0, meanR, twZMat, tZMat, TZvar, LZvar, fitResultF, fitResultX, fitResultFz, fitResultXz, tDMatrix, alphas_a, tg_frame] = analyzeShapeDynamics(cellInfo, tInt, T0, pixSize, divededVarEachFrm, PathNameResults, diamThresh, smthThresh)
% Function to calculate, plot, fit data from measureShapeDynamics
% Calculates:
% assymetry/constriction site position
% curvature of the cell backbone
% Fits models to constriction data
%
% -----
% Input:
% cellInfo: structure with data from measureShapeDynamics (different cells and
% frames)
%
% tInt: interval time between successive frames in minutes
%
% T0: time after synchrony for the first frame.
%
% pixSize: pixel size in nm (30nm/px for N-SIM)
%
% divededVarEachFrm: matrix with info whether cells are divided from measureShapeDynamics
%
% PathNameResults: Path where to store results, string example 'C:\folder\folder'
%
% -------
% Output:
% twMatrix: Matrix containing width versus time info for all cells:
% [t wDiam] (relative width)
%
% tlMatrix: Matrix containing length versus time info for all cells: [t length]
%
% Tvar: Matrix with Tc and Tg for each cell (measured by different models)
%
% Lvar: Matrix with Lb, Lc, LG for each cell (measured by different models)
%
% BIC: Matrix with Bayesian information criterion: statistical test to
% asses whether the Xiao model is overfitting.
%
% p1C0: Matrix containing the position of the constriction site, measured
% by the MTS channel, for each time point for reach cell.
%
% meanR: Matrix of mean radius of curvature for each timepoint for each
% cell.
%
% twZMat: Matrix containing the diameter at the constriction site measured
% by MTS, the width of the (primary) Z-ring or spot and the time, for each
% time point for each cell: [t, diameter by MTS, diameter by Z]
%
% tZMat:Matrix containing information on the state of the Zring for each
% timepoint for each cell: [t, position of Z, relative Intensity at midcell, Zassembled?]
% 
% TZvar: Similar to Tvar, but with Tc fixed in the model fits.
% TZvar = [TcFz' TgFz' TcXz' TgXz'];
% 
% LZvar: Similar to Lvar, but with Tc fixed in the model fits.
% LZvar = [LcFz' LgFz' Lb' LcXz' LgXz'];
% 
% fitResultF, fitResultX, fitResultFz, fitResultXz: variables containing
% the fit results (fitted parameters and summed normalized residuals) for
% each cell (corresponding to the rows).
%  
% tDMatrix: Matrix containing waist diameter and maximal diameter vs time
% for each cell. For N cells:
% frame 1:   [t1 DiamMTS1 DiamMTSmax1 t2 DiamMTS2 DiamMTSmax2 ... tN DiamMTSN DiamMTSmaxN]
% frame 2:   [t1 DiamMTS1 DiamMTSmax1 t2 DiamMTS2 DiamMTSmax2 ... tN DiamMTSN DiamMTSmaxN]
% frame ...: [...                     ...                     ... ...                    ]
% last frame:[t1 DiamMTS1 DiamMTSmax1 t2 DiamMTS2 DiamMTSmax2 ... tN DiamMTSN DiamMTSmaxN]
% 
% alphas_a: Fitted alpha's form Xiao model with onset time fixed to Tc.
% 
% tg_frame: First frame in which the cell is divided.

%% init param

[frNum, ~]= size(cellInfo);
[~, cellNum] = size(cellInfo{1,1});
% improve speed: preallocation of matrices
twMatrix(frNum,cellNum*2)=0;
twZMat(frNum,cellNum*3)=0;
tZMat(frNum,cellNum*4)=0;
tlMatrix(frNum,cellNum*2)=0;
c0Matrix(frNum,cellNum*2)=0;
p1C0(frNum,cellNum)=0;
p2C0(frNum,cellNum)=0;
meanR(frNum,cellNum)=0;
tForR(frNum,cellNum)=0;
lForR(frNum,cellNum)=0;
Tz(1,cellNum)=0;
ZintAll{cellNum,frNum}=[];
stdOvermean=zeros(frNum,cellNum);
tDMatrix(frNum,cellNum*3)=0;

try %Check if dual colour data.
    cellInfo{1, 1}{1, 1}.Zassembled;
    DualC=true;
catch
    DualC=false;
end

%% Toggles which plots to plot etc.
togglePlots

%% init Curvature msrmnt figure
if plotCurv
    curvFig=figure;
    title('Radius of curvature against normalized time, all cells.')
end

%% 1st loop through cells and frames: calculate/create matrices
% twMatrix  = time and normalized waist width for each cell+frame
%       cell1  cell2 ... celln
% fr1   t1 w1  t2 w2 ... tn wn
% fr2
% tlMatrix  = time and length for each cell+frame
%       cell1  cell2 ... celln
% fr1   t1 l1  t2 l2 ... tn ln
% fr2
% c0Matrix  = constriction position x,y for each cell+frame
%       cell1       cell2      ... celln
% fr1   c0x1 c0y1   c0x2 c0y2  ... c0xn c0yn
% fr2
% twZMat    = time, absolute diameter by MTS and FtsZ resp. for each cell+frame
%       cell1           cell2           ... celln
% fr1   t1 dMTS1 dFtsZ1 t2 dMTS2 dFtsZ2 ... tn dMTSn dFtsZn
% fr2
% tZMat     = time, relative position of the Zring wrt the long axis of the
% cell, relative intensity in the Z channel at the midcell region,
% Zassembled? (1=yes, 0=no)
%       cell1                             ... celln
% fr1   t1 Zpos1 midcellInt1 Zassembled1  ... tn Zposn midcellIntn Zassembledn
% fr2
% p1C0 & p2C0 = relative position of constriction (min MTS diam) along cell
% axis. p1C0>p2C0, p1C0+p2C0=1.
%       cell1   cell2   ... celln
% fr1   t1 pos1 t2 pos2 ... tn posn
% fr2
% Tz        = vector containing the time FtsZ assembled for each cell

for cellIdx = 1 : cellNum
    % Calculate curvature, constriction position, construct twMatrix, tDMatrix, tlMatrix, tZMat, ZintAll, c0Matrix
    for frIdx = 1 :  frNum
        
        divededVarTemp = divededVarEachFrm(divededVarEachFrm(:,1) == frIdx, 2:3);
        
        % If exists or is divided
        if (~isempty(cellInfo{frIdx, 1}{1, cellIdx}) && isfield(cellInfo{frIdx, 1}{1, cellIdx},'c0')) || divededVarTemp(cellIdx, 2) == 1
            
            % If it is not divided
            if divededVarTemp(cellIdx, 2) == 0
                
                % Grab the position c0 along the bacteria where the diam
                % is minimum
                c0 = cellInfo{frIdx, 1}{1, cellIdx}.c0;
                % Grab the waist diam
                wDiam = cellInfo{frIdx, 1}{1, cellIdx}.waistD;
                DiamMTS=cellInfo{frIdx, 1}{1, cellIdx}.minD.*pixSize;
                DiamMTSmax=cellInfo{frIdx, 1}{1, cellIdx}.maxD.*pixSize;
                if DualC
                    try
                        diamZ=cellInfo{frIdx, 1}{1, cellIdx}.diamZ;
                        Zpos=cellInfo{frIdx, 1}{1, cellIdx}.Zpos;
                        midcellInt=cellInfo{frIdx, 1}{1, cellIdx}.midcellInt;
                        Zassembled=cellInfo{frIdx, 1}{1, cellIdx}.Zassembled;
                    catch
                        diamZ=NaN;
                        Zpos=NaN;
                        midcellInt=NaN;
                        Zassembled=false;
                    end
                end
                % Grab the length of the bacteria
                leng = cellInfo{frIdx, 1}{1, cellIdx}.bactLength*pixSize;
                t = (frIdx-1)*tInt+T0;
                % Grab mesh of the bacteria
                mesh = cellInfo{frIdx, 1}{1, cellIdx}.mesh;
                
                twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t wDiam];
                tDMatrix(frIdx, cellIdx*3 - 2 : cellIdx*3 ) = [t DiamMTS DiamMTSmax];
                tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t leng];
                if DualC
                    twZMat(frIdx, cellIdx*3 - 2 : cellIdx*3 ) = [t DiamMTS diamZ];
                    tZMat(frIdx, cellIdx*4 - 3 : cellIdx*4 ) = [t Zpos midcellInt 0];
                    if Zassembled
                        tZMat(frIdx, cellIdx*4)=1;
                    end
                end
                
                if calcMTSmid || calcCurv % Create (smoothed) centerline
                    % Centerline, pole position and curvature stuff
                    % Find the bacteria central line
                    cLine = zeros(size(mesh,1),2);
                    cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
                    cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);
                    
                    nBreak = 10;
                    % smooth the centreLine
                    [xC, yC] = spline2d(cLine(:,1),cLine(:,2),nBreak);
                end
                
                if calcMTSmid % Calculate constriction position by MTS channel
                    c0Matrix(frIdx, cellIdx*2 - 1 : cellIdx*2) = c0;
                    % Find the point of the central line closest to the center c0Cell
                    distancesFromC = sqrt((xC - c0(1, 1)).^2 + (yC - c0(1, 2)).^2);
                    [~, indexOfMin] = min(distancesFromC);
                    
                    % Compute bacteria pole-center length normalize by the tot lenght
                    totLength = arclength( xC, yC, 'spline');
                    p1CenterLength = arclength( xC(1 : indexOfMin), yC(1 : indexOfMin), 'spline')/totLength;
                    p2CenterLength = arclength( xC(indexOfMin : end), yC(indexOfMin : end), 'spline')/totLength;
                    
                    % Save p1C and p2C accordingly with their length
                    if p1CenterLength >= p2CenterLength
                        p1C0(frIdx, cellIdx) = p1CenterLength;
                        p2C0(frIdx, cellIdx) = p2CenterLength;
                    else
                        p1C0(frIdx, cellIdx) = p2CenterLength;
                        p2C0(frIdx, cellIdx) = p1CenterLength;
                    end
                end
                
                if calcCurv % Calculate curvature of centerline
                    % Compute the curvature (R from lsq fit)for each cell
                    if plotCurv
                        figure(curvFig),
                    end
                    try
                        meanR(frIdx, cellIdx) = computeCurv(cellInfo{frIdx, 1}{1, cellIdx}.c0, cellInfo{frIdx, 1}{1, cellIdx}.mesh, pixSize, plotCurv);
                    catch curvError
                        throw(curvError)
                    end
                    tForR(frIdx, cellIdx) = t;
                    lForR(frIdx, cellIdx) = leng;
                    
                    if plotCurv
                        % Debug:
                        figure(curvFig),
                        subplot(1, 2, 2)
                        %       imagesc(img),
                        plot(mesh(:,1), mesh(:,2),'y-', 'LineWidth', 2);hold on
                        plot(mesh(:,3), mesh(:,4),'y-', 'LineWidth', 2);
                        plot(xC, yC, 'b-')
                        plot(c0(1,1),c0(1,2),'r*')
                        plot(xC(1), yC(1), 'ro')
                        plot(xC(end), yC(end), 'ro')
                        %                       legend('meshLeft', 'meshRight', 'centralLine', 'centerDiv', 'pole1', 'pole2' )
                        axis equal
                        hold off
                        drawnow
                    end
                end
                
                % ZintAll
                try
                    ZintAll{cellIdx,frIdx}=cellInfo{frIdx, 1}{1, cellIdx}.ZintLength;
                catch
                end
                
            else % if it is divided
                
                t = (frIdx-1)*tInt+T0;
                
                twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t 0];
                tDMatrix(frIdx, cellIdx*3 - 2 : cellIdx*3 ) = [t 0 0];
                if DualC
                    twZMat(frIdx, cellIdx*3 - 2 : cellIdx*3 ) = [t 0 0];
                    tZMat(frIdx, cellIdx*4 - 3 : cellIdx*4 ) = [t 0 0 0];
                end
                tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [t nan];
                c0Matrix(frIdx, cellIdx*3 - 1 : cellIdx*3) = [nan nan];
            end
        else
            twMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [nan nan];
            tDMatrix(frIdx, cellIdx*3 - 2 : cellIdx*3 ) = [nan nan nan];
            tlMatrix(frIdx, cellIdx*2 - 1 : cellIdx*2 ) = [nan nan];
            c0Matrix(frIdx, cellIdx*2 - 1 : cellIdx*2) = [nan nan];
        end
        
    end % end frame cycle
    
    if DualC % determine when Z first assembles
        %             if ~FtsW
        Tz(cellIdx)=getTz(tZMat(:,cellIdx*4),tZMat(:,(cellIdx*4)-3));
        %             else
        %                 % do the midcell thing and get the std/mean and arrival y/n
        %                 tTemp=T0:tInt:(T0+(frNum-1)*tInt);
        %                 for frIdx=1:frNum
        %                     try
        %                         midCellStart=0.3*length(ZintAll{cellIdx,frIdx});
        %                         midCellStop=0.7*length(ZintAll{cellIdx,frIdx});
        %                         midCellInt=ZintAll{cellIdx,frIdx}(:,round(midCellStart):round(midCellStop));
        %                         stdOvermean(frIdx,cellIdx)=std(midCellInt(2,:))/mean(midCellInt(2,:));
        %                     catch
        %                     end
        %                 end
        %                 Tz(cellIdx)=getTz(stdOvermean(:,cellIdx)>0.29,tTemp);
        %             end
    end
    
    % If the cell is present only in less than 5 frames
    if length( find( ( isnan(twMatrix(:, cellIdx*2 - 1 )) ) == 0 )) < 5 || sum( find( twMatrix(:, cellIdx*2)  > 0 ) ) < 5
        twMatrix(:, cellIdx*2 - 1 : cellIdx*2 ) = repmat([nan nan], size(twMatrix, 1), 1);
        tlMatrix(:, cellIdx*2 - 1 : cellIdx*2 ) = repmat([nan nan], size(twMatrix, 1), 1);
        c0Matrix(:, cellIdx*2 - 1 : cellIdx*2) = repmat([nan nan], size(twMatrix, 1), 1);
    end
end% end cell cycle



% Initialize var
wNoZeroTransf = [];
tNoZeroTransf = [];
wNoZeroTransfF = [];
tNoZeroTransfF = [];
lNoZeroTransf = [];
lNoZeroTransfF = [];
tlNoZeroTransf = [];
tlNoZeroTransfF = [];
normLNoZeroTransf = [];
mkdir([PathNameResults '\matlabPlot\LvsTCell']);
wNoZeroTransfX=[];
tNoZeroTransfX=[];
wNoZeroTransfXa=[];
tNoZeroTransfXa=[];
wNoZeroTransfX24=[];
tNoZeroTransfX24=[];

% initialize figures
if plotdWvT
    deltaW=figure;
    title('constriction aka delta D vs time')
end
if DualC && plotSCWvTZ
    allZposNorm=figure;
    hold on, grid on
    title('Relative position of the Z-ring, all cells.')
    if plotCurvAll
        allMidMTS=figure;
        hold on, grid on
        title('Relative position of the minimal MTS diameter, all cells.')
    end
end
if plotCurvAll && calcCurv
    curvAll=figure;
    hold on, grid on
    xlabel('Normalized Time (t-tc)/(tg-tc)')
    ylabel('Radius of curvature[nm]')
    legend('Curvature as a function of time')
    legend('Location', 'best')
end


% Improve speed: preallocation
Lc2D(cellNum) = 0;
Lg2D(cellNum) = 0;
LcF(cellNum) = 0;
LgF(cellNum) = 0;
Lb(cellNum) = 0;
Tc2D(cellNum) = 0;
TcF(cellNum) = 0;
Tg2D(cellNum) = 0;
TgF(cellNum) = 0;
LcX(cellNum) = 0;
LgX(cellNum) = 0;
TcX(cellNum) = 0;
TgX(cellNum) = 0;
Lz(cellNum) = 0;
tg_frame(cellNum)= 0;
alphas_a(cellNum)= 0;
TcSmth(cellNum)  = 0;
TcSpline(cellNum)= 0;
LcSmth(cellNum)  = 0;
LcSpline(cellNum)= 0;

if fitSCZ_Feingold
    TcFz(cellNum) = 0;
    TgFz(cellNum) = 0;
    LcFz(cellNum) = 0;
    LgFz(cellNum) = 0;
end
if fitSCZ_Xiaofree
    TcXz(cellNum) = 0;
    TgXz(cellNum) = 0;
    LcXz(cellNum) = 0;
    LgXz(cellNum) = 0;
end

BICF=zeros(1,cellNum);
BICX=zeros(1,cellNum);
% Good cells are the one present for more than 10 frames
cellIdxGood = 1;
fitResultF=zeros(4,cellNum);
fitTcWidth=zeros(2,cellNum);
fitResultX=zeros(5,cellNum);
fitResultFz=zeros(4,cellNum);
fitResultXz=zeros(5,cellNum);

%% 2nd loop through cells and frames: plot figures with fits
%
for cellIdx = 1 : cellNum
    
    % Measure Lc adn Lg: length at Tc when constriction start and Tg when the
    % division occurs
    
    % Delete each raw with a w = nan
    twMatrixTemp = twMatrix(:, 2*cellIdx-1 : 2*cellIdx);
    twMatrixTemp = twMatrixTemp(~any(isnan(twMatrixTemp),2),:);
    twMatrixTempWithZero = twMatrixTemp(~any(isnan(twMatrixTemp),2),:);
    % Delete each raw with w = 0
    twMatrixTemp(any(twMatrixTemp == 0, 2), :) = [];
    
    % Delete each raw with a l = nan
    tlMatrixTemp = tlMatrix(:, 2*cellIdx-1 : 2*cellIdx);
    tlMatrixTemp = tlMatrixTemp(~any(isnan(tlMatrixTemp),2),:);
    % Delete each raw with l = 0
    tlMatrixTemp(any(twMatrixTemp == 0, 2), :) = [];
    
    % If the cell has any information
    if length(twMatrixTemp(:, 1)) > 1
        
        % No zero
        tTemp = twMatrixTemp(:, 1);
        wTemp =  twMatrixTemp(:, 2);
        tlTemp = tlMatrixTemp(:, 1);
        lTemp =  tlMatrixTemp(:, 2);
        
        % With zero
        tTempWithZero = twMatrixTempWithZero(:, 1);
        wTempWithZero = twMatrixTempWithZero(:, 2);
        
        lengTemp =  tlMatrixTemp(:, 2);
        
        % FIT w vs time with 2D and feingold model
        % with zero, to include the point where the cell is divided for
        % sure.
        [tc, tg, ~, fitOutTimeNoBinNoLeng, ~, res2D]  = fitSeamusModBound(tTempWithZero, wTempWithZero);
        [tcF, tgF, wMaxF, fitOutTimeNoBinNoLengF, resnormF, resF, BICf]  = fitFeingoldModBound(tTempWithZero, wTempWithZero);
        
        % Find first time for which Dmax-Dmin > threshold
        tc_id=find(tDMatrix(:, 3*cellIdx)-tDMatrix(:, 3*cellIdx-1) > diamThresh ,1);
        tc_d=tDMatrix( tc_id , 3*cellIdx-2); %
        if ~isempty(tc_d)
            lc_d = mean(lengTemp(tc_d-tInt <= tTemp & tTemp <= tc_d+tInt));
        else
            tc_d=NaN;
            lc_d=NaN;
        end
        
        % Fit model to obtain length function
        [tcX, tgX, wMaxX, alpha, fitOutTimeNoBinNoLengX, resnormX, resX, BICx]  = fitXiao(tTempWithZero, wTempWithZero); % fit with Jie Xiao model
        
        %% Smoothing data for threshold
        wSmth=smooth(wTempWithZero,3); % 3 seems good
        ppwSpline=splinefit(tTempWithZero,wTempWithZero,20); %ppX = splinefit(t,x,nBreak); % 20 breakpoints is better than 10
        tTempUpsampled=linspace(min(tTempWithZero),max(tTempWithZero),500);
        wSpline=ppval(ppwSpline,tTempUpsampled);
        TcSmth(cellIdx)=tTempWithZero(find(wSmth>max(wSmth)*smthThresh,1,'last'));
        TcSpline(cellIdx)=tTempUpsampled(find(wSpline>max(wSpline)*smthThresh,1,'last'));
        LcSmth(cellIdx) = mean(lengTemp(TcSmth(cellIdx)-tInt <= tTemp & tTemp <= TcSmth(cellIdx)+tInt));
        LcSpline(cellIdx)= mean(lengTemp(TcSpline(cellIdx)-tInt <= tTemp & tTemp <= TcSpline(cellIdx)+tInt));
        tTempSmth=tTempWithZero;
        
        %%
        if fitSC_Xiaofree_Tc_Tg || (fitSC_XiaoAlpha && ~isnan(tc_d))
            tg_id=find(twMatrixTemp(:,2),1,'last');
            if tg_id>tc_id+4
                tg_frame(cellIdx)= twMatrixTemp( find(twMatrixTemp(:,2),1,'last') ,1) +tInt./2;
                wTempWithZeroFilt=wTempWithZero(tc_id:tg_id);
                tTempWithZeroFilt=tTempWithZero(tc_id:tg_id);
%                 "Resolution" by removing low diameter datapoints. Causes
%                 trouble by not leaving enough points for fitting
%                 tTempWithZeroFilt=tTempWithZeroFilt(wTempWithZeroFilt>0.5);
%                 wTempWithZeroFilt=wTempWithZeroFilt(wTempWithZeroFilt>0.5);
            elseif tg_id<4
                wTempWithZeroFilt=wTempWithZero;
                tTempWithZeroFilt=tTempWithZero;
            else
                wTempWithZeroFilt=wTempWithZero(1:tg_id);
                tTempWithZeroFilt=tTempWithZero(1:tg_id);
            end
            try
                [tcXa, tgXa, wMaxXa, alpha_a, fitOutTimeNoBinNoLengXa, resnormXa, resXa, ~]  = fitXiao(tTempWithZeroFilt, wTempWithZeroFilt, 'Tc', TcSpline(cellIdx)); % Constrain Tc and Tg %,'Tg', tg_frame(cellIdx)
            catch
                disp('ddd')
            end
            alphas_a(cellIdx)=alpha_a;
            
        end
        if DualC
            if fitSC_Feingold_Tz
                [tcFz, tgFz, wMaxFz, fitOutTimeNoBinNoLengFz, resnormFz, ~, ~]  = fitFeingoldModBound(tTempWithZero, wTempWithZero,Tz(cellIdx));
                fitResultFz(:,cellIdx)=[tcFz, tgFz, wMaxFz, resnormFz];
            else
                fitResultFz=[];
            end
            if fitSC_Xiaofree_Tz
                [tcXz, tgXz, wMaxXz, alphaz, fitOutTimeNoBinNoLengXz, resnormXz, ~, ~]  = fitXiao(tTempWithZero, wTempWithZero,'Tc',TcSpline(cellIdx)); % Tz(cellIdx)
                fitResultXz(:,cellIdx)=[tcXz, tgXz, wMaxXz, alphaz, resnormXz];
            elseif fitSC_Xiaofree_Tc_Tg
                [tcXz, tgXz, wMaxXz, alphaz, fitOutTimeNoBinNoLengXz, resnormXz, ~, ~]  = fitXiao(tTempWithZero, wTempWithZero,'Tc',TcSpline(cellIdx),'Tg',tg_frame(cellIdx));
                fitResultXz(:,cellIdx)=[tcXz, tgXz, wMaxXz, alphaz, resnormXz];
            else
                fitResultXz=[];
            end
        end
        
        BICF(cellIdx)=BICf;
        BICX(cellIdx)=BICx;
        
        fitResultF(:,cellIdx)=[tcF, tgF, wMaxF, resnormF];
        fitResultX(:,cellIdx)=[tcX, tgX, wMaxX, alpha, resnormX];
        fitTcWidth(:,cellIdx)=[tc_d lc_d];
        
        idxCellWithATimeIntAroundTc = find(tc-tInt <= tTemp & tTemp <= tc+tInt);
        idxCellWithATimeIntAroundTg = find(tg - 3*tInt <= tTemp & tTemp <= tg);
        idxCellWithATimeIntAroundTcF = find(tcF-tInt <= tTemp & tTemp <= tcF+tInt);
        idxCellWithATimeIntAroundTgF = find(tgF - 3*tInt <= tTemp & tTemp <= tgF);
        idxCellWithATimeIntAroundTcX = find(tcX-tInt <= tTemp & tTemp <= tcX+tInt);
        idxCellWithATimeIntAroundTgX = find(tgX - 3*tInt <= tTemp & tTemp <= tgX);
        idxOfCellBetweenTcTg = find(tc <= tTemp & tTemp <= tg);
        
        if DualC
                idxCellWithATimeIntAroundTcFz = find(Tz(cellIdxGood) - tInt <= tTemp & tTemp <= Tz(cellIdxGood)+tInt);
                Lz(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTcFz));
        end
        
        Lc2D(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTc));
        Lg2D(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTg));
        LcF(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTcF));
        LgF(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTgF));
        LcX(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTcX));
        LgX(cellIdxGood) = mean(lengTemp(idxCellWithATimeIntAroundTgX));
        Tc2D(cellIdxGood) = tc;
        TcF(cellIdxGood) = tcF;
        Tg2D(cellIdxGood) = tg;
        TgF(cellIdxGood) = tgF;
        TcX(cellIdxGood) = tcX;
        TgX(cellIdxGood) = tgX;
        
        % FIT l vs time, model: L = l0.*exp(kL.*t)
        [kL, l0, fitOutLength]  = fitLengthVsTime(tlTemp, lTemp);
        
        %         % FIT log(l) vs time, model: L = log(l0log) + kL.*t;
        %         [kLlog, l0log, fitOutLengthlog]  = fitLengthVsTimeLog(tlTemp, log(lTemp));
        %
        %         % FIT normalize l vs time, model: L = nl0.*exp(kL.*t)
        %         [nkL, nl0, nfitOutLength]  = fitNormLengthVsTime(tlTemp, lTemp./lTemp(end));
        
        % Use elongation fit to guess Lg, Lc and Lb: lOut = l0.*exp(kL.*t);
        %         Lc2D(cellIdxGood) = l0.*exp(kL.*tc);
        %         Lg2D(cellIdxGood) = l0.*exp(kL.*tg);
        %         LcF(cellIdxGood) = l0.*exp(kL.*tcF);
        %         LgF(cellIdxGood) = l0.*exp(kL.*tgF);
        %         LcX(cellIdxGood) = l0.*exp(kL.*tcX);
        %         LgX(cellIdxGood) = l0.*exp(kL.*tgX);
        Lb(cellIdxGood) = l0;
        
        % Transformation of t => (t-tc)/(tg-tc)
        if cellIdx == 1
            tWithZeroTransf=zeros(frNum,1);
            wWithZeroTransf=zeros(frNum,1);
        end
        
        try
            tWithZeroTransf(:, cellIdxGood) = (tTempWithZero - tc)/(tg -tc);
            %tlWithZeroTransf(:, cellIdx) = (tlTempWithZero - tc)/(tg -tc);
            wWithZeroTransf(:, cellIdxGood) = wTempWithZero;
            %lWithZeroTransf(:, cellIdx) = lTempWithZero;
        catch
            disp('Matrix dimensions do not fit, let me fix it for you!')
            missing = size(tWithZeroTransf,1)- size(tTempWithZero,1);
            if missing > 0
                fix = linspace((frNum-1)*tInt+T0,(frNum-1)*tInt+T0,missing);
                fix=rot90(fix);
                tTempWithZero=[tTempWithZero ; fix];
                wTempWithZero=[wTempWithZero;fix*0];
            else
                fprintf('missing is %d!\n',missing)
                
            end
            try
                tWithZeroTransf(:, cellIdx) = (tTempWithZero - tc)/(tg -tc);
                wWithZeroTransf(:, cellIdx) = wTempWithZero;
            catch
                disp('Damaged beyond repair, debug time!')
                % output for debug
                tTempWithZero %#ok<NOPRT>
                size(tTempWithZero)
                wTempWithZero %#ok<NOPRT>
                size(wTempWithZero)
                size(tc)
                size(tg)
            end
        end
        
        lengthWtemp = length(wTemp);
        lNoZeroTransf(end + 1: end + lengthWtemp, 1) = lTemp;
        normLNoZeroTransf(end + 1: end + lengthWtemp, 1) = lTemp./lTemp(end);
        wNoZeroTransf(end + 1: end + lengthWtemp, 1) = wTemp;
        tNoZeroTransf(end + 1: end + lengthWtemp, 1) = (tTemp - tc)/(tg - tc);
        tlNoZeroTransf(end + 1: end + lengthWtemp, 1) = (tlTemp - tc)/(tg - tc);
        
        % Transformation of t => (t-tcF)/(tgF-tcF)
        tWithZeroTransfF(:, cellIdxGood) = (tTempWithZero - tcF)/(tgF -tcF);
        %tlWithZeroTransfF(:, cellIdx) = (tlTempWithZero - tcF)/(tgF -tcF);
        wWithZeroTransfF(:, cellIdxGood) = wTempWithZero;
        %lWithZeroTransfF(:, cellIdx) = lTempWithZero;
        
        %         tWithZeroTransfX(:, cellIdxGood) = (tTempWithZero - tcX)/(tgX -tcX);
        %         wWithZeroTransfX(:, cellIdxGood) = wTempWithZero;
        
        lengthWtemp = length(wTemp);
        lNoZeroTransfF(end + 1: end + lengthWtemp, 1) = lTemp;
        wNoZeroTransfF(end + 1: end + lengthWtemp, 1) = wTemp;
        tNoZeroTransfF(end + 1: end + lengthWtemp, 1) = (tTemp - tcF)/(tgF - tcF);
        tlNoZeroTransfF(end + 1: end + lengthWtemp, 1) = (tlTemp - tcF)/(tgF - tcF);
        
        %%
        
        if fitSC_Xiaofree
            wNoZeroTransfX(end + 1: end + lengthWtemp, 1) = wTemp;
            tNoZeroTransfX(end + 1: end + lengthWtemp, 1) = (tTemp - tcX)/(tgX - tcX);
        end
        if fitSC_XiaoAlpha
            wNoZeroTransfXa(end + 1: end + lengthWtemp, 1) = wTemp;
            tNoZeroTransfXa(end + 1: end + lengthWtemp, 1) = (tTemp - tcX)/(tgX - tcX);
        end
        
        % plot delta W vs time
        if plotdWvT
            dWtemp=wTemp(2:end)-wTemp(1:end-1);
            figure(deltaW),
            hold on,
            plot((tTemp(2:end) - tcF)./(tgF - tcF),dWtemp,'o')
        end
        
        % Plot l vs time for each cell with its own fitting
        if plotSCLvT
            figure
            plot(fitOutLength(:, 1), fitOutLength(:, 2), 'r')
            hold on
            % plot(tTemp, wTemp, 'bo')
            plot(tlTemp, lTemp, 'bo')
            plot(0,l0,'+')
            
            text(60, 2500,['$$ F_{model}(t) = '  num2str(l0,'%.0f') '\exp{( ' num2str(kL) ' *t) } $$'],'interpreter','latex', 'FontSize',12 )
            text(60, 2800, '$$ F_{model}(t) = l_{0}*\exp{( k_{L}*t )}$$' ,'interpreter','latex', 'FontSize',12 )
            set(gca,'box','on')
            
            title(['Length vs time. Cell ' num2str(cellIdx)], 'FontSize', 20, 'FontName', 'Calibri Light'),
            legend ('Fitting with exponential model', 'Length', 'Binned data', 'Location', 'best')
            xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
            ylabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light' );
            hold off
            
            saveas(gcf, [PathNameResults '\matlabPlot\LvsTCell' num2str(cellIdx) '.fig'])
        end
        
        if calcMTSmid
            % Take the distance pole-center of division for the cell with a
            % diameter w < 0.8
            p1C0GoodCell = p1C0(idxOfCellBetweenTcTg, cellIdx);
            p2C0GoodCell = p2C0(idxOfCellBetweenTcTg, cellIdx);
            
            %         p1C0GoodCellAll{cellIdx} =  p1C0GoodCell;
            %         p2C0GoodCellAll{cellIdx} =  p2C0GoodCell;
            if plotCentPos
                figure,
                plot(tTemp(idxOfCellBetweenTcTg),  p1C0GoodCell, 'r*'), hold on,
                plot(tTemp(idxOfCellBetweenTcTg),  p2C0GoodCell, 'b*'),
                xlabel('Time [min]')
                ylabel('Pole-Center Distance normalize by the tot lenght')
                legend('Pole1-Center Distance', 'Pole2-Center Distance')
                legend('Location', 'best')
            end
        end
        % Plot the curvature as a function of normalized time (t-tc)/(tg-tc)
        if plotCurv
            figure,
            plot((tForR(meanR(:,cellIdx)~=0, cellIdx)-tc)/(tg-tc), meanR(meanR(:,cellIdx)~=0, cellIdx), 'k')
            xlabel('Normalized Time (t-tc)/(tg-tc)')
            ylabel('Curvature (Radius of the Osculating Circle)[nm]')
            legend('Curvature as a function of time')
            legend('Location', 'best')
        end
        
        % Plot w vs time for each cell with its own fitting
        if plotSCWvT
            figure
            hold on
            plot(tTempWithZero, wTempWithZero, 'bo:')
            plotSCWvTLeg{1,1}='Data';
            if fitSC_Feingold
                plot(fitOutTimeNoBinNoLengF(:, 1), fitOutTimeNoBinNoLengF(:, 2), 'r', 'DisplayName', 'Feingold fit')
                % put formula on graph
                text(50, 0.8,['$$ F_{model}(t) = '  num2str(wMaxF,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tcF)) ' }{' num2str(round(tgF)) ' - '  num2str(round(tcF)) '} \right )^2}$$'],'interpreter','latex', 'Color', 'r' ,'FontSize',12 )
            end
            if fitSC_Xiaofree
                plot(fitOutTimeNoBinNoLengX(:, 1), fitOutTimeNoBinNoLengX(:, 2),'g', 'DisplayName', 'Xiao fit')
                % put formula on graph
                text(50, 0.9, ['$$ F_{model}(t) = ' num2str(wMaxX,'%.2f') '({1 -  \left ( \frac{ t - ' num2str(round(tcX)) ' }{' num2str(round(tgX)) ' - '  num2str(round(tcX)) '})^{' num2str(alpha) '} \right ) )^{1/' num2str(alpha) '}}$$'],'interpreter','latex','Color','g', 'FontSize',12 )
            end
            if fitSC_Xiaofree_Tc_Tg %tcXz, tgXz, wMaxXz, alphaz, resnormXz
                plot(fitOutTimeNoBinNoLeng(:, 1), XiaoModel(fitOutTimeNoBinNoLeng(:, 1),tcXz, tgXz, wMaxXz, alphaz) ,'g--', 'DisplayName', 'Xiao fit, fixed Tc and Tg')
                text(50, 0.9, ['$$ F_{model}(t) = ' num2str(wMaxXz,'%.2f') '({1 -  \left ( \frac{ t - ' num2str(round(tcXz)) ' }{' num2str(round(tgXz)) ' - '  num2str(round(tcXz)) '})^{' num2str(alphaz) '} \right ) )^{1/' num2str(alphaz) '}}$$'],'interpreter','latex','Color','g', 'FontSize',12 )
            elseif fitSC_XiaoAlpha
%                 plot(fitOutTimeNoBinNoLengXa(:, 1), fitOutTimeNoBinNoLengXa(:, 2),'g--', 'DisplayName', 'Xiao fit, fixed Tc and Tg')
                plot(fitOutTimeNoBinNoLeng(:, 1), XiaoModel(fitOutTimeNoBinNoLeng(:, 1),tcXa,tgXa,wMaxXa,alpha_a) ,'g--', 'DisplayName', 'Xiao fit, fixed Tc and Tg')
                % put formula on graph
                text(50, 0.9, ['$$ F_{model}(t) = ' num2str(wMaxXa,'%.2f') '({1 -  \left ( \frac{ t - ' num2str(round(tcXa)) ' }{' num2str(round(tgXa)) ' - '  num2str(round(tcXa)) '})^{' num2str(alpha_a) '} \right ) )^{1/' num2str(alpha_a) '}}$$'],'interpreter','latex','Color','g', 'FontSize',12 )
            end
            if fitSC_2D
                plot(fitOutTimeNoBinNoLeng(:, 1), fitOutTimeNoBinNoLeng(:, 2),'b', 'DisplayName', 'Fit 2D model')
            end
            if DualC
                if fitSC_Feingold_Tz
                    plot(fitOutTimeNoBinNoLengFz(:, 1), fitOutTimeNoBinNoLengFz(:, 2), 'm', 'DisplayName', 'Feingold fit with fixed Tc')
                    % put formula on graph
                    text(50, 0.6,['$$ F_{model}(t) = '  num2str(wMaxFz,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tcFz)) ' }{' num2str(round(tgFz)) ' - '  num2str(round(tcFz)) '} \right )^2}$$'],'interpreter','latex', 'Color', 'm' ,'FontSize',12 )
                end
                if fitSC_Xiaofree_Tz
                    plot(fitOutTimeNoBinNoLengXz(:, 1), fitOutTimeNoBinNoLengXz(:, 2),'c', 'DisplayName', 'Fit Xiao, fixed Tc')
                    % put formula on graph
                    text(50, 0.5, ['$$ F_{model}(t) = ' num2str(wMaxXz,'%.2f') '({1 -  \left ( \frac{ t - ' num2str(round(tcXz)) ' }{' num2str(round(tgXz)) ' - '  num2str(round(tcXz)) '})^{' num2str(alphaz) '} \right ) )^{1/' num2str(alphaz) '}}$$'],'interpreter','latex','Color','c', 'FontSize',12 )
                end
            end
            if plotSC_fitRes
                if fitSC_Feingold
                    plot(fitOutTimeNoBinNoLengF(:, 1), resF,'r', 'DisplayName', 'residuals Feingold fit')
                end
                if fitSC_Xiaofree
                    plot(fitOutTimeNoBinNoLengX(:, 1), resX,'g', 'DisplayName', 'residuals Xiao fit')
                end
                if fitSC_2D
                    plot(fitOutTimeNoBinNoLeng(:, 1), res2D,'b', 'DisplayName', 'residuals 2D model')
                end
            end
            
            plot(tTempSmth, wSmth, 'b-', 'DisplayName', 'Sliding average')
            plot(tTempUpsampled, wSpline, 'k-', 'DisplayName', 'Spline')
            plot([0 300],[smthThresh*max(wSpline), smthThresh*max(wSpline)],'k--', 'DisplayName', 'Threshold on spline')
            title(['Normalized waist width vs time. Cell ' num2str(cellIdx)], 'FontSize', 20, 'FontName', 'Calibri Light'),
            legend('show')
            legend('Location', 'best')
            set(gca,'box','on')
            xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
            ylabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );
            grid on
            hold off
            
            saveas(gcf, [PathNameResults '\matlabPlot\WvsTCell' num2str(cellIdx) '.fig'])
        end
        
        if DualC && plotSCWvTZ
            figure,
            hold on
            plot(twZMat(:,cellIdx*3 - 2),twZMat(:,cellIdx*3-1), 'DisplayName', 'MTS diameter') %  DiamMTS vs t
            plot(twZMat(:,cellIdx*3 - 2),twZMat(:,cellIdx*3), 'DisplayName', 'FtsZ diameter') % diamZ vs t
            if fitSCZ_Feingold
                if sum(tZMat(:, cellIdx*4))>4 % Check Zassembled==1 at enough points or not.
                    % fit diamZ vs t with Feingold model, but only points where
                    % The peak for the Z channel is at midcell
                    try
                        [tcFZ, tgFZ, ~, fitOutTimeNoBinNoLengFZ, ~,  ~, ~]=fitFeingoldModBound...
                            (twZMat((tZMat(:, cellIdx*4-2)>0.3).*(tZMat(:, cellIdx*4-2)<0.7)==1 ,cellIdx*3 - 2) , ...
                            twZMat( (tZMat(:, cellIdx*4-2)>0.3).*(tZMat(:, cellIdx*4-2)<0.7)==1 ,cellIdx*3)); % This is getting ridiculously long
                    catch
                        tcFZ=NaN;
                        tgFZ=NaN;
                        fitOutTimeNoBinNoLengFZ=[NaN NaN];
                    end
                    plot(fitOutTimeNoBinNoLengFZ(:,1),fitOutTimeNoBinNoLengFZ(:,2), 'DisplayName','fit Feingold')
                    TcFz(cellIdxGood)=tcFZ;
                    TgFz(cellIdxGood)=tgFZ;
                    LcFz(cellIdxGood)=l0.*exp(kL.*tcFZ);
                    LgFz(cellIdxGood)=l0.*exp(kL.*tgFZ);
                end
            end
            if fitSCZ_Xiaofree
                if sum(tZMat(:, cellIdx*4))>4
                    % fit diamZ vs t with Xiao model, but only points where
                    % The peak for the Z channel is at midcell
                    try
                        [tcXZ, tgXZ, ~, ~, fitOutTimeNoBinNoLengXZ, ~,  ~, ~]=fitXiao...
                            (twZMat((tZMat(:, cellIdx*4-2)>0.3).*(tZMat(:, cellIdx*4-2)<0.7)==1 ,cellIdx*3 - 2) , ...
                            twZMat( (tZMat(:, cellIdx*4-2)>0.3).*(tZMat(:, cellIdx*4-2)<0.7)==1 ,cellIdx*3));
                    catch
                        tcXZ=NaN;
                        tgXZ=NaN;
                        fitOutTimeNoBinNoLengXZ=[NaN NaN];
                    end
                    plot(fitOutTimeNoBinNoLengXZ(:,1),fitOutTimeNoBinNoLengXZ(:,2), 'DisplayName','fit Xiao')
                      TcXz(cellIdxGood)=tcXZ;
                    TgXz(cellIdxGood)=tgXZ;
                    LcXz(cellIdxGood)=l0.*exp(kL.*tcXZ);
                    LgXz(cellIdxGood)=l0.*exp(kL.*tgXZ);
                end
            end
            if fitSCmts_Feingold
                % fit diameter MTS vs t with Feingold model
                [tcFZ, tgFZ, wMaxFZ, fitOutTimeNoBinNoLengFZ, ~,  ~, ~]=fitFeingoldModBound(twZMat(:,cellIdx*3 - 2) , twZMat(:,cellIdx*3-1));
                plot(fitOutTimeNoBinNoLengFZ(:,1),fitOutTimeNoBinNoLengFZ(:,2), 'DisplayName', 'fit Feingold MTS')
%                 legendsSCWvTZ=[legendsSCWvTZ,];
                TcFz(cellIdxGood)=tcFZ;
                TgFz(cellIdxGood)=tgFZ;
                LcFz(cellIdxGood)=l0.*exp(kL.*tcFZ);
                LgFz(cellIdxGood)=l0.*exp(kL.*tgFZ);
            end
            if fitSCmts_Xiaofree
                % fit diameter MTS vs t with Xiao model
                [tcXZ, tgXZ, wMaxXZ, alphaXZ, fitOutTimeNoBinNoLengXZ, ~,  ~, ~]=fitXiao(twZMat(:,cellIdx*3 - 2) , twZMat(:, cellIdx*3-1));
                plot(fitOutTimeNoBinNoLengXZ(:,1),fitOutTimeNoBinNoLengXZ(:,2), 'DisplayName','fit Xiao MTS')
                TcXz(cellIdxGood)=tcXZ;
                TgXz(cellIdxGood)=tgXZ;
                LcXz(cellIdxGood)=l0.*exp(kL.*tcXZ);
                LgXz(cellIdxGood)=l0.*exp(kL.*tgXZ);
            end
            
            title(['Diameter vs time for cell ' num2str(cellIdx)])
            legend('Location', 'best')
            saveas(gcf, [PathNameResults '\matlabPlot\FtsZ_WvsTCell' num2str(cellIdx) '.fig'])
            % Let's check Zpos and flip if needed
            ZposTest=tZMat(tZMat(:, cellIdx*4-2)~=0, cellIdx*4-2);
            % flip it if...
            % if FtsZ was below 0.2 and never above 0.8
            % if it never was above 0.8 or below 0.2 AND ends below 0.5 and didnt start
            if ( ~any(ZposTest>0.8) && any(ZposTest<0.2) ) || ( ~any(ZposTest>0.8) && ~any(ZposTest<0.2) && ZposTest(end)<0.5 )
                ZposTest=1-ZposTest;
            end
            figure(allZposNorm), hold on,
            plot( ( tZMat(tZMat(:,cellIdx*4-2)~=0,cellIdx*4 - 3) -tc)/(tg-tc) ,ZposTest) %Zpos normalized time (no 0's)
            if plotCurvAll
                figure(allMidMTS), hold on,
                plot( ( tZMat(tZMat(:,cellIdx*4-2)~=0,cellIdx*4 - 3) -tc)/(tg-tc) ,p1C0(tZMat(:,cellIdx*4-2)~=0,cellIdx)) % MTS mid normalized time (no 0's)
            end
            figure,
            hold on
            grid on
            plot(tZMat(tZMat(:, cellIdx*4-2)~=0,cellIdx*4 - 3),ZposTest)%Zpos
            plot(tZMat(:,cellIdx*4 - 3),tZMat(:, cellIdx*4-1))
            plot(tZMat(:,cellIdx*4 - 3),tZMat(:, cellIdx*4))
            plot(tZMat(:,cellIdx*4 - 3),p1C0(:,cellIdx))
            title(['Z ring intensity, Z ring position and invagination position for cell ' num2str(cellIdx)])
            legend('position of Z (max)', 'Z midcell Intensity', 'Z assembled?', 'Constriction position by MTS')
            legend('Location', 'best')
            saveas(gcf, [PathNameResults '\matlabPlot\ZringvsTCell' num2str(cellIdx) '.fig'])
            
            if plotCurvAll && calcCurv %plot curvature for all cells vs normalized time
                figure(curvAll),
                plot((tForR(meanR(:,cellIdx)~=0, cellIdx)-tc)/(tg-tc), meanR(meanR(:,cellIdx)~=0, cellIdx))
            end
        end
        
        cellIdxGood = cellIdxGood +1;
    else
        disp('The cell')
        disp(cellIdx)
        disp('has been deleted')
    end
    
end % end cycle over each cell
if plotCurvAll
    figure(allMidMTS)
    saveas(gcf, [PathNameResults '\assymMTS.fig'])
end
if DualC && plotSCWvTZ
    figure(allZposNorm)
    saveas(gcf, [PathNameResults '\assymFtsZ.fig'])
end

%%

if exist('BICX13','var') && exist('BICX24','var')
    BIC=[BICF;BICX;BICX13;BICX24];
else
    BIC=[BICF;BICX];
end

% For each cell save Tc , Tg, etc
% if not calculated, make empty variables to not mess up Tvar and Lvar
% assignment
if ~fitSC_2D
    Tc2D=zeros(1,cellNum);
    Tg2D=zeros(1,cellNum);
    Lc2D=zeros(1,cellNum);
    Lg2D=zeros(1,cellNum);
end
if ~fitSC_Feingold
    TcF=zeros(1,cellNum);
    TgF=zeros(1,cellNum);
    LcF=zeros(1,cellNum);
    LgF=zeros(1,cellNum);
end
if ~fitSC_Xiaofree
    TcX=zeros(1,cellNum);
    TgX=zeros(1,cellNum);
    LcX=zeros(1,cellNum);
    LgX=zeros(1,cellNum);
end
Tvar = [TcF' TgF' TcX' TgX' Tc2D' Tg2D' Tz' fitTcWidth(1,:)' TcSmth' TcSpline'];
Lvar = [LcF' LgF' Lb' LcX' LgX' Lc2D' Lg2D' Lz' fitTcWidth(2,:)' LcSmth' LcSpline'];

if ~fitSCZ_Feingold
    TcFz=zeros(1,cellNum);
    TgFz=zeros(1,cellNum);
    LcFz=zeros(1,cellNum);
    LgFz=zeros(1,cellNum);
end
if ~fitSCZ_Xiaofree
    TcXz=zeros(1,cellNum);
    TgXz=zeros(1,cellNum);
    LcXz=zeros(1,cellNum);
    LgXz=zeros(1,cellNum);
end
TZvar = [TcFz' TgFz' TcXz' TgXz'];
LZvar = [LcFz' LgFz' Lb' LcXz' LgXz'];

%% -----------------------------------------------------------------------%
% PLOT SECTION All CELL TOGETHER                                          %
%-------------------------------------------------------------------------%

%     figure,
%     plot(Lg - Lc, Lc, 'r*'), hold on
%     xlabel('Lc - Lg - length at the division - length at the costriction [nm]')
%     ylabel('Lc - length at the costriction [nm]')

% 2D area modell all cell with also cell already divided--------------
if fitSC_2D
    tTrasf = reshape(tWithZeroTransf, [], 1);
    %tlTrasf = reshape(tlWithZeroTransf, [], 1);
    wTrasf = reshape(wWithZeroTransf, [], 1);
    %lTrasf = reshape(lWithZeroTransf, [], 1);
    
    [tc, tg, wMax, fitOutTimeNoBinNoLeng]  = fitSeamusModBound( tTrasf, wTrasf);
    
    figure
    plot(fitOutTimeNoBinNoLeng(:, 1), fitOutTimeNoBinNoLeng(:, 2), 'r')
    hold on
    plot(tTrasf, wTrasf, 'bo')
    
    text(1.5, 0.8,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )}$$'],'interpreter','latex', 'FontSize',12 )
    text(1.5, 0.9, '$$ F_{model}(t) = w*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )}$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title('Caulobacter septum waist ratio r/R vs time. All Cell.', 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with 2D septum area model all cells', 'Caulobacter during division')
    xlabel('Time [min]','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',18, 'FontName', 'Calibri Light' );
end
%% Feingold modell all cell with also cell already divided
if fitSC_Feingold
    tTrasfF = reshape(tWithZeroTransfF, [], 1);
    %tlTrasf = reshape(tlWithZeroTransf, [], 1);
    wTrasfF = reshape(wWithZeroTransfF, [], 1);
    %lTrasf = reshape(lWithZeroTransf, [], 1);
    [tcF, tgF, wMaxF, fitOutTimeNoBinNoLengF]  = fitFeingoldModBound( tTrasfF, wTrasfF); % [tc, tg, wMax, fitOut, resnorm, residual, BIC]
    
    figure
    plot(fitOutTimeNoBinNoLengF(:, 1), fitOutTimeNoBinNoLengF(:, 2), 'r')
    hold on
    plot(tTrasfF, wTrasfF, 'bo')
    
    text(1.5, 0.8,['$$ F_{model}(t) = '  num2str(wMaxF,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tcF)) ' }{' num2str(round(tgF)) ' - '  num2str(round(tcF)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
    text(1.5, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title('Caulobacter septum waist ratio r/R vs time. All Cell.', 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with Feingold model all cells', 'Caulobacter during division', 'Residuals')
    xlabel('Time [min]','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',18, 'FontName', 'Calibri Light' );
    grid on
end
%% %Waist vs time for all cells with real time
%% parameter renaming (Why? Ask Anna)

tTrasfNoZero = tNoZeroTransf ;
tlTrasfNoZero = tlNoZeroTransf ;
%wTrasfNoZero = reshape(wNoZeroTransf, [], 1);
wTrasfNoZero = wNoZeroTransf;
lTrasfNoZero = lNoZeroTransf;
normLTrasfNoZero = normLNoZeroTransf;

if fitSC_Feingold
    tTrasfNoZeroF = tNoZeroTransfF ;
    tlTrasfNoZeroF = tlNoZeroTransfF ;
    %wTrasfNoZero = reshape(wNoZeroTransf, [], 1);
    wTrasfNoZeroF = wNoZeroTransfF;
    lTrasfNoZeroF = lNoZeroTransfF;
end
if fitSC_Xiaofree
    tTrasfNoZeroX = tNoZeroTransfX;
end
%% 2D model normalized

% Delete negative value %also needed for AllLvT
negTimeIdx = find( tTrasfNoZero<0);
tTrasfNoZeroNoNeg = tTrasfNoZero;
tTrasfNoZeroNoNeg(negTimeIdx) = [];

if fitSC_2D
    wTrasfNoZeroNoNeg = wNoZeroTransf;
    wTrasfNoZeroNoNeg(negTimeIdx) = [];
    
    normLTrasfNoZeroNoNeg = normLNoZeroTransf;
    normLTrasfNoZeroNoNeg(negTimeIdx) = [];
    
    [tc, tg, wMax, fitOutTimeNoBinNoLeng, ~, res2D]  = fitSeamusModBound( tTrasfNoZeroNoNeg, wTrasfNoZeroNoNeg);
    
    figure
    plot(fitOutTimeNoBinNoLeng(:, 1), fitOutTimeNoBinNoLeng(:, 2), 'r')
    hold on
    plot( tTrasfNoZeroNoNeg ,  wTrasfNoZeroNoNeg, 'bo')
    plot(tTrasfNoZeroNoNeg, res2D,'o')
    [resTsort2D, I]=sort(tTrasfNoZeroNoNeg);
    resSrt2D=res2D(I);
    lag=ceil(length(resSrt2D)/20);
    resavg2D=tsmovavg(resSrt2D,'s',lag,1);
    plot(resTsort2D,resavg2D)
    
    text(0.5, 0.8,['$$ F_{model}(t) = '  num2str(wMax,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tc)) ' }{' num2str(round(tg)) ' - '  num2str(round(tc)) '} \right )}$$'],'interpreter','latex', 'FontSize',12 )
    text(0.5, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )}$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title('Waist width vs normalized time. All Cells.', 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with 2D septum area model all cells', 'Data with rescaled time', 'Residuals')
    xlabel('Time [min]','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',18, 'FontName', 'Calibri Light' );
    grid on
end
%% FIT l vs time
if plotSCLvT && plotAllLvT
    lTrasfNoZeroNoNeg = lNoZeroTransf;
    lTrasfNoZeroNoNeg(negTimeIdx) = [];
    [kL, l0, fitOutLength, ~, resL]  = fitLengthVsTime( tTrasfNoZeroNoNeg, lTrasfNoZeroNoNeg);
    % Plot l vs time for each cell with its own fitting
    figure
    plot(fitOutLength(:, 1), fitOutLength(:, 2), 'r')
    hold on
    plot(tTrasfNoZeroNoNeg ,  lTrasfNoZeroNoNeg, 'bo')
    plot(fitOutLength(:, 1), resL)
    
    text(0.6, 3000,['$$ F_{model}(t) = '  num2str(l0,'%.0f') '\exp{ ' num2str(kL) ' *t} $$'],'interpreter','latex', 'FontSize',12 )
    text(0.6, 2600, '$$ F_{model}(t) = l_{0}*\exp{ k_{L}*t }$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title( 'Caulobacter length vs time. Population ' , 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with exponential model', 'Length data', 'Residuals')
    xlabel('Time [min]','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Length [nm]','FontSize',18, 'FontName', 'Calibri Light' );
    grid on
    hold off
end
%% FIT normalize l vs time
if plotSCLvT && plotAllLvT && plotAllnormLvT
    [nkL, nl0, nfitOutLength]  = fitNormLengthVsTime( tTrasfNoZeroNoNeg, normLTrasfNoZeroNoNeg);
    % Plot l vs time for each cell with its own fitting
    figure
    plot(nfitOutLength(:, 1), nfitOutLength(:, 2), 'r')
    hold on
    plot(tTrasfNoZeroNoNeg ,  normLTrasfNoZeroNoNeg, 'bo')
    
    text(0.6, 0.9,['$$ F_{model}(t) = '  num2str(nl0,'%.0f') '\exp{ ' num2str(nkL) ' *t} $$'],'interpreter','latex', 'FontSize',12 )
    text(0.6, 0.6, '$$ F_{model}(t) = l_{0}*\exp{ k_{L}*t }$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title( 'Length vs time. All cells.' , 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with exponential model', 'Caulobacter during division', 'Binned data')
    xlabel('Time [min]','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Length [nm]','FontSize',18, 'FontName', 'Calibri Light' );
    hold off
end
%% Feingold model normalized
if fitSC_Feingold
    % Delete negative value
    negTimeIdxF = find( tTrasfNoZeroF<0);
    tTrasfNoZeroNoNegF = tTrasfNoZeroF;
    tTrasfNoZeroNoNegF(negTimeIdxF) = [];
    
    wTrasfNoZeroNoNegF = wNoZeroTransfF;
    wTrasfNoZeroNoNegF(negTimeIdxF) = [];
    
    [tcF, tgF, wMaxF, fitOutTimeNoBinNoLengF, ~, resF, BICf]  = fitFeingoldModBound( tTrasfNoZeroNoNegF, wTrasfNoZeroNoNegF);
    
    figure
    plot(fitOutTimeNoBinNoLengF(:, 1), fitOutTimeNoBinNoLengF(:, 2), 'r')
    hold on
    plot( tTrasfNoZeroNoNegF ,  wTrasfNoZeroNoNegF, 'bo')
    plot(tTrasfNoZeroNoNegF, resF,'o')
    [resTsortF, I]=sort(tTrasfNoZeroNoNegF);
    resSrtF=resF(I);
    lag=50;
    if lag >= resSrtF
        lag=1;
    end
    resavgF=tsmovavg(resSrtF,'s',lag,1);
    plot(resTsortF,resavgF)
    
    text(0.5, 0.8,['$$ F_{model}(t) = '  num2str(wMaxF,'%.2f') '\sqrt{1 -  \left ( \frac{ t - ' num2str(round(tcF)) ' }{' num2str(round(tgF)) ' - '  num2str(round(tcF)) '} \right )^2}$$'],'interpreter','latex', 'FontSize',12 )
    text(0.5, 0.9, '$$ F_{model}(t) = wMax*\sqrt{1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^2}$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title('Waist width vs normalized time. All Cells.', 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with Feingold model all cells', 'Data with rescaled time', 'Residuals')
    xlabel('Normalized time, (t-tc)/(tg-tc)','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Normalized waist width','FontSize',18, 'FontName', 'Calibri Light' );
    grid on
end
%% Xiao model free fit normalized and fit again with Xiao
if fitSC_Xiaofree
    % Delete negative value
    negTimeIdxX = find( tTrasfNoZeroX<0);
    tTrasfNoZeroNoNegX = tTrasfNoZeroX;
    tTrasfNoZeroNoNegX(negTimeIdxX) = [];
    
    wTrasfNoZeroNoNegX = wNoZeroTransfX;
    wTrasfNoZeroNoNegX(negTimeIdxX) = [];
    
    [tcX, tgX, wMaxX, alpha, fitOutTimeNoBinNoLengX, ~, resX,BICx]  = fitXiao(tTrasfNoZeroNoNegF, wTrasfNoZeroNoNegF); % Normalize with Feingold or Xiao
    
    figure
    plot(fitOutTimeNoBinNoLengX(:, 1), fitOutTimeNoBinNoLengX(:, 2), 'r')
    hold on
    plot(tTrasfNoZeroNoNegF,wTrasfNoZeroNoNegF, 'bo')
    plot(tTrasfNoZeroNoNegF, resX,'o')
    % moving average?
    [resTsortX, I]=sort(tTrasfNoZeroNoNegF);
    resSrtX=resX(I);
    lag=50;
    if lag >= resSrtX
        lag=1;
    end
    resavgX=tsmovavg(resSrtX,'s',lag,1);
    plot(resTsortX,resavgX)
    
    text(0.3, 0.9, ['$$ F_{model}(t) = ' num2str(wMaxX,'%.2f') '({1 -  \left ( \frac{ t - ' num2str(round(tcX)) ' }{' num2str(round(tgX)) ' - '  num2str(round(tcX)) '})^{' num2str(alpha) '} \right ) )^{1/' num2str(alpha) '}}$$'],'interpreter','latex', 'FontSize',12 )
    text(0.5, 0.9, '$$ F_{model}(t) = w_{Max}*(1 - \left (  \frac{t - t_c }{ t_g - t_c} \right )^{alpha})^{1/alpha}$$' ,'interpreter','latex', 'FontSize',12 )
    set(gca,'box','on')
    
    title('Caulobacter septum waist ratio r/R vs normalized time for Xiao fit with free Alpha. All Cells', 'FontSize', 23, 'FontName', 'Calibri Light'),
    legend ('Fitting with Xiao model', 'Data with rescaled time', 'Residuals')
    xlabel('Time normalized with Feingold tc and tg','FontSize',18, 'FontName', 'Calibri Light'  );
    ylabel('Septum diameter/max diameter(along bacteria length)','FontSize',18, 'FontName', 'Calibri Light' );
    grid on
end
%% histograms of BIC and alpha
if plotBIChist
    figure,
    hist(BICF)
    title('Histogram BIC Feingold')
    meanBICF=mean(BICF);
    stdBICF=std(BICF);
    text(meanBICF,1,['Mean= ' num2str(meanBICF) ', standard deviation = ' num2str(stdBICF)])
    figure,
    hist(BICX)
    title('Histogram BIC Xiao')
    meanBICX=mean(BICX);
    stdBICX=std(BICX);
    text(meanBICX,1,['Mean= ' num2str(meanBICX) ', standard deviation = ' num2str(stdBICX)])
end

figure,hist(fitResultX(4,:))
title('Histogram of alphas from fits to single cells')
meanAlpha=mean(fitResultX(4,:));
stdAlpha=std(fitResultX(4,:));
text(meanAlpha,1,['Mean= ' num2str(meanAlpha) ', standard deviation = ' num2str(stdAlpha)])

if DualC && fitSC_Xiaofree_Tz
    figure,hist(fitResultXz(4,:))
    title('Histogram of alphas from fits to single cells, with constrained Tc')
    meanAlpha=mean(fitResultXz(4,:));
    stdAlpha=std(fitResultXz(4,:));
    text(meanAlpha,1,['Mean= ' num2str(meanAlpha) ', standard deviation = ' num2str(stdAlpha)])
end
end
%------------------------------------------------------------------------%
% OTHER METHODS
%------------------------------------------------------------------------%

function [xOut,yOut]=spline2d(x,y,nBreak)
nPt = numel(x);
nPtOut = 10*nPt;% upsample the output, but keep discrete

t=linspace(0,1,nPt);
tOut=linspace(0,1,nPtOut);
ppX = splinefit(t,x,nBreak);
ppY = splinefit(t,y,nBreak);
xOut = ppval(ppX,tOut);
yOut = ppval(ppY,tOut);
end
%------------------------------------------------------------------------%

function [R] = computeCurv(c0Cell, mesh, pixSize, plotCurv)

% Grab the central line
cLine = zeros(size(mesh,1),2);
cLine(:,1) = mean([mesh(:,1),mesh(:,3)],2);
cLine(:,2) = mean([mesh(:,2),mesh(:,4)],2);

nBreak = 10;
% smooth the centreLine
[xCLine, yCLine] = spline2d(cLine(:,1),cLine(:,2),nBreak);
%cLineSmooth=[xC(:),yC(:)];

[k, dfdx, d2fdx2, x, y] = getCurvature(xCLine, yCLine);

% Grab the mean curvature
goodCurvIdx = find(k < mean(k) + std(k) & k > mean(k) - std(k));
meanK = mean(k(goodCurvIdx));

% Find the point of the central line closest to the center c0Cell
distancesFromC = sqrt((xCLine - c0Cell(1, 1)).^2 + (yCLine - c0Cell(1, 2)).^2);
[minDistance, indexOfMin] = min(distancesFromC);

xCirc = (0:pi/64:2*pi);

% Plot the central line with the osculating circle
if plotCurv
    subplot(1, 2, 1)
    plot(mesh(:,1), mesh(:,2),'y-', 'LineWidth', 2);hold on
    plot(mesh(:,3), mesh(:,4),'y-', 'LineWidth', 2);
    plot(xCLine, yCLine, 'b-')
    plot(xCLine(indexOfMin), yCLine(indexOfMin), 'r*')
    plot(c0Cell(1,1), c0Cell(1,2), 'c*')
    plot(x, y', 'g+')
    xlabel('[nm]')
    ylabel('[nm]')
end
% Plot osculating circle
Rk = 1/meanK;
d2fdx2Mean =  mean(d2fdx2);
dfdxMean = mean(dfdx);
if d2fdx2Mean > 0
    if  dfdxMean > 0
        Xc = xCLine(indexOfMin) - Rk*sin( atan(dfdxMean) );
        Yc = yCLine(indexOfMin) + Rk*cos( atan(dfdxMean) );
        %                     plot(Rk*cos(xCirc) + Xc, Rk*sin(xCirc) + Yc, 'r--')
        %                     axis equal
    else
        Xc = xCLine(indexOfMin) - Rk*sin( atan(dfdxMean) );
        Yc = yCLine(indexOfMin) + Rk*cos( atan(dfdxMean) );
        %                     plot(Rk*cos(xCirc) + Xc, Rk*sin(xCirc) + Yc, 'r--')
        %                     axis equal
    end
    
else
    if dfdxMean > 0
        Xc = xCLine(indexOfMin) + Rk*sin( atan(dfdxMean) );
        Yc = yCLine(indexOfMin) - Rk*cos( atan(dfdxMean) );
        %                     plot(Rk*cos(xCirc) + Xc, Rk*sin(xCirc) + Yc, 'r--')
        %                     axis equal
    else
        Xc = xCLine(indexOfMin) + Rk*sin( atan(dfdxMean) );
        Yc = yCLine(indexOfMin) - Rk*cos( atan(dfdxMean) );
        %                     plot(Rk*cos(xCirc) + Xc, Rk*sin(xCirc) + Yc, 'r--')
        %                     axis equal
        
    end
end
%             legend('meshLeft', 'meshRight', 'centralLine', 'Division Center', 'Division Center','centralLineRotInsideFunc', 'Osculating circle' )
%             axis equal
%             hold off

xCguess = Xc;
yCguess = Yc;
Rguess = Rk;
[Rtemp, xC, yC]  = fitCurveWithCircle(xCLine, yCLine, xCguess, yCguess, Rguess);
R = Rtemp*pixSize;

% Plot circle from fit
if plotCurv
    plot(Rtemp*cos(xCirc) + xC, Rtemp*sin(xCirc) + yC, 'r--')
end
%             legend('meshLeft', 'meshRight', 'centralLine', 'Division Center', 'Division Center','centralLineRotInsideFunc', 'Osculating circle' )
axis equal
hold off
end
%------------------------------------------------------------------------%

function wOut = XiaoModel(t,tc,tg,wMax,alpha)
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
