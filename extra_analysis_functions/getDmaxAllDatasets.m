for analysisCycleId=1:17 % play with this
    %% ------------
    tInt=5;
    pixSize=30;
    switch analysisCycleId
        case 1 %ZWTA
            T0=22;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160622 DC WT t0 22min\Results\together 3-5_9_10_12-15\testRerun2';
        case 2 %ZWTB
            T0=30;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160713 DC WT t0 30min\Results\together_05-12\testRerun2';
        case 3 %ZMutA
            T0=80;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160613 DC ML2159 t0 80min\Results\together 1_2_4-9\testRerun';
        case 4 %ZMutB
            T0=35;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160621 DC FtsZ ML159 35min\together 2_4-11_13_15\testRerun';
        case 5 %WWTA
            T0=50;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160816 WT DC 1305 FtsWGFP t0 50min\Results\together 01-06 and 10\manual Tw 2';
        case 6 %WWTB
            T0=37;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160926 WT FtsW 18 488 37 min\Results\Together 01-07\manual Tw 2';
        case 7 %WMutA
            T0=36; 
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160927 DC ML2159 FtsW t0 36min\Results\together_01-10\manual Tw 3';
        case 8 %WMutC
            T0=34;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\161005 FtsW Mut 34min\Results\together 1-12\manual Tw 2';
        case 9 % ZfWTA FOM
            T0=80; tInt=10;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160526 dc fOSFO T080\Results\together 1-9, 11,12,14-17';
        case 10 % NWTA
            T0=45; tInt=5;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160727 WT FTsN 45min\Results\together 1-11';
        case 11 %NWTB
            T0=28;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170412 WT DC FtsN t28\Results\together 1-4 11-15';
        case 12 %NMutA
            T0=50;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160726 ML2159 DC FtsN t0 50min\Results\together 1_6_7_9';
        case 13 %NMutB
            T0=25;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\160810 ML2159 FtsN t0 25 min\Results\ML2159 FtsN t25 01_06_07_09 together';
        case 14 % ZpWT
            T0=22;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170221 WT PYE T0 22min\Results\together 1-5_7_9_11\Cleaned';
        case 15 % PYE Mut
            T0=16;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170307 ML2159 DC FtsZ PYE\Results\together 1-3_5_6_8_14_15';
        case 16 % ZfWTB FOM
            T0=17; tInt=10;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170412 WT DC FtsZ fosfo t17\Results\together 1-5_12-15';
        case 17 % ZfMutA FOM
            T0=25; tInt=10;
            PathNameResults='\\lebsrv1.epfl.ch\LEB\SHARED\Bacteria Shape Project\SIM data\Analysis\170704 ML2159 FosfoT0 25 min\Results\together_1-7';
    end
    %%
    
    load([PathNameResults '\cellInfoWithC0.mat'])
    load([PathNameResults '\testRerun_Thr_100_sThr_0.92\Tvar.mat']);
    DmaxAll=getDmaxAll(cellInfo);
    save([PathNameResults '\DmaxAll.mat'],'DmaxAll')
    tTemp=T0:tInt:T0+74*tInt;
    Tc=Tvar(:,10);
    for cellIdx=size(DmaxAll,1):-1:1
        DmaxAtTcAll(cellIdx,1)=nanmean(DmaxAll(cellIdx,Tc(cellIdx,1)-tInt <= tTemp & tTemp <= Tc(cellIdx,1)+tInt),2);
    end
    save([PathNameResults '\DmaxAtTcAll.mat'],'DmaxAtTcAll')
    clear DmaxAll DmaxAtTcAll cellInfo
end