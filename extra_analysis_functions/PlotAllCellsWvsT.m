% Plot w vs time for each cell on the same figure
num_frames=size(cellInfo,1);
num_cells=size(cellInfo{1,1},2);
% close all
T=0:tInt:(num_frames-1)*tInt; % initialize time vector
T=T+T0;
% initialize figures
WvT = figure; 
WvL = figure;
LvT = figure;
WvLvT= figure;
% Events=figure;
WvTfits=figure;
LvTfits=figure;
WvLvTfits=figure;

% ASAtest=figure;

% preallocate matrices
WMat=zeros(num_cells,num_frames);
LMat=zeros(num_cells,num_frames);
Whist=zeros(num_cells*num_frames,2);
Lhist=zeros(num_cells*num_frames,2);
ElMat=zeros(num_cells,2);
for cellIdx=1:num_cells;
    W=zeros(1,num_frames); % preallocate vectors
    L=zeros(1,num_frames);
    for frIdx=1:num_frames;
        try
            W(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD;
            WMat(cellIdx,frIdx)=cellInfo{frIdx,1}{1,cellIdx}.minD/cellInfo{frIdx,1}{1,cellIdx}.maxD; 
            Whist(frIdx+(cellIdx-1)*num_frames,:)=[(frIdx-1)*5,W(frIdx)];
            L(frIdx)=cellInfo{frIdx,1}{1,cellIdx}.bactLength*30;
            LMat(cellIdx,frIdx)=cellInfo{frIdx,1}{1,cellIdx}.bactLength*30;
            Lhist(frIdx+(cellIdx-1)*num_frames,:)=[L(frIdx),W(frIdx)];
        catch
            W(frIdx)=0;
            WMat(cellIdx,frIdx)=0;
            Whist(frIdx+(cellIdx-1)*num_frames,:)=[(frIdx-1)*tInt, NaN];
            try
                if W(frIdx-1)==0;
                    L(frIdx)=NaN;
                    LMat(cellIdx,frIdx)=NaN;
                    Lhist(frIdx+(cellIdx-1)*num_frames,:)=[NaN NaN];
                else
                    L(frIdx)=L(frIdx-1);
                    LMat(cellIdx,frIdx)=L(frIdx-1);
                    Lhist(frIdx+(cellIdx-1)*num_frames,:)=[L(frIdx-1), 0];
                end
            catch
                L(frIdx)=NaN;
                LMat(cellIdx,frIdx)=NaN;
                Lhist(frIdx+(cellIdx-1)*num_frames,:)=[NaN NaN];
            end
        end
%         figure(ASAtest),
%         hold on
%         try
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),cellInfo{frIdx,1}{1,cellIdx}.diffSum,'+r') % diffSum
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),cellInfo{frIdx,1}{1,cellIdx}.areaPrev,'+g') % areaPrev
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),cellInfo{frIdx,1}{1,cellIdx}.areaNow,'+b') % areaNow
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),cellInfo{frIdx,1}{1,cellIdx}.diffSum/cellInfo{frIdx,1}{1,cellIdx}.areaPrev,'+k') %diffSum norm
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),cellInfo{frIdx,1}{1,cellIdx}.areaNow-cellInfo{frIdx,1}{1,cellIdx}.areaPrev,'+c') % delta area
%             plot(size(cellInfo{frIdx,1}{1,cellIdx}.recordSteps),(cellInfo{frIdx,1}{1,cellIdx}.areaNow-cellInfo{frIdx,1}{1,cellIdx}.areaPrev)/cellInfo{frIdx,1}{1,cellIdx}.areaPrev,'+m') % delta area norm
%         catch
% %             fprintf('Failed for cell %d in frame %d\n', cellIdx, frIdx)
%         end
    end
    
    
    figure(WvT),
    hold on
    plot(T,W,'.:','MarkerSize',15); %,'.','MarkerSize',15
    figure(WvL),
    hold on
    plot(L,W,'-','LineWidth', 0.5); % ,'.','MarkerSize',15
    figure(LvT),
    hold on
    plot(T,L,'-','LineWidth', 0.5);
%     figure(WvLvT), % plot trajectories and key events
%     hold on
%     plot3(T,L,W,'-','LineWidth', 0.5)
%     plot3(0,Lb(cellIdx),1,'+','MarkerSize',10) % birth
%     plot3(Tc(cellIdx),Lc(cellIdx),1,'o','MarkerSize',10) % constriction onset
%     plot3(Tg(cellIdx),Lg(cellIdx),0,'.','MarkerSize',20) % division
%     figure(Events), % plot only key events
%     hold on
%     plot3(0,Lvar(cellIdx,5),1,'+','MarkerSize',10) % birth
%     plot3(Tvar(cellIdx,3),Lvar(cellIdx,3),1,'o','MarkerSize',10) % constriction onset
%     plot3(Tvar(cellIdx,4),Lvar(cellIdx,4),0,'.','MarkerSize',20) % division
%     plot3([0, Tvar(cellIdx,3), Tvar(cellIdx,4)],[Lvar(cellIdx,5), Lvar(cellIdx,3), Lvar(cellIdx,4)],[1, 1, 0], ':') % plot a line to connect them
    
    figure(WvTfits), % WvT fitting and plotting
    hold on
    [tcF, tgF, wMaxF, fitOutTimeNoBinNoLengF]  = fitFeingoldModBound(T, W);
    plot(fitOutTimeNoBinNoLengF(:, 1), fitOutTimeNoBinNoLengF(:, 2), '-')
    figure(LvTfits), % LvT fitting and plotting
    hold on
    Lf=L;
    Tf=T;
    for i=size(L,2):-1:1
        if isnan(Lf(i))
            Lf(i)=[];
            Tf(i)=[];
        end
    end
    [kL, l0, fitOutLength]  = fitLengthVsTime(Tf, Lf);
    ElMat(cellIdx,:)=[kL, l0];
    plot(fitOutLength(:,1), fitOutLength(:,2))
%     figure(WvLvTfits),
%     hold on
%     plot3(fitOutLength(:,1),fitOutLength(:,2),fitOutTimeNoBinNoLengF(1:size(fitOutLength,1),2))
%     plot3(0,Lb(cellIdx),1,'+','MarkerSize',10) % birth
%     plot3(Tc(cellIdx),Lc(cellIdx),1,'o','MarkerSize',10) % constriction onset
%     plot3(Tg(cellIdx),Lg(cellIdx),0,'.','MarkerSize',20) % division
end

figure(WvT),
title('Normalized waist width vs time. All Cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

figure(WvL),
title('Normalized waist width vs Length. All Cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

figure(LvT),
title('Length vs time. All Cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light' );

figure(WvLvT),
title('Normalized waist width vs Length vs Time. All Cells', 'FontSize', 20, 'FontName', 'Calibri Light'),
xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
ylabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light'  );
zlabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

% figure(Events)
% title('Birth, constriction onset and division', 'FontSize', 20, 'FontName', 'Calibri Light'),
% xlabel('Time [min]','FontSize',16, 'FontName', 'Calibri Light'  );
% ylabel('Length [nm]','FontSize',16, 'FontName', 'Calibri Light'  );
% zlabel('Normalized waist width','FontSize',16, 'FontName', 'Calibri Light' );

% Calculate and plot the mean width for each timepoint
Wmean=sum(WMat); % There was an issue with the mean function.
Wmean=Wmean./num_cells;
figure(WvT),
plot(T,Wmean,'k','LineWidth', 3)
% Lmean=sum(LMat); % doesn't work with NaN...
% Lmean=Lmean./num_cells;
% figure(LvT),
% plot(T,Lmean,'k','LineWidth', 3)

% % Plot a 3D histogram of W vs time
% figure,
% hist3(Whist,[20,50])
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
       
LhistOrd=sortrows(Lhist);
% figure, plot(LhistOrd(:,1),LhistOrd(:,2))
Lsmth=smooth(LhistOrd(:,1),LhistOrd(:,2),100);
figure(WvL), plot(LhistOrd(:,1),Lsmth,'k-', 'LineWidth', 2)

figure, plot(ElMat(:,2),ElMat(:,1),'o')
title('Elongation rate vs Lb', 'FontSize', 20, 'FontName', 'Calibri Light'),


 
        