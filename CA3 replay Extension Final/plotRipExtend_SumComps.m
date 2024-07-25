%%% plot ripple extension summary plots
% Plots comparison of shuffles for the ripple extension model
% LKW 9/5/2021
set(0,'DefaultLineLineWidth',2)

distNoiseFlag = 1;
opsinNoiseFlag = 1;
VNoiseFlag = 0;

pStruct.nRamps = 21;
pStruct.nPs = 101;
    
%% Light Scattering Noise Comparisons
if distNoiseFlag == 1    
    % Load Variables
    % No Distance Noise
    FRC_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IMA_1000_CI95.mat');
    FRP_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IP_1000_CI95.mat');
    DRC_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IMA_1000_CI95.mat');
    DRP_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IP_1000_CI95.mat');
    BRC_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IMA_1000_CI95.mat');
    BRP_NDN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IP_1000_CI95.mat');
    % Distance Noise, 8mW Source, 5mW Threshold
    FRC_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\FR_IMA_DNoise_Irr8_perms1000_CI95.mat');
    FRP_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\FR_IP_DNoise_Irr8_perms1000_CI95.mat');
    DRC_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\DR_IMA_DNoise_Irr8_perms1000_CI95.mat');
    DRP_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\DR_IP_DNoise_Irr8_perms1000_CI95.mat');
    BRC_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\BR_IMA_DNoise_Irr8_perms1000_CI95.mat');
    BRP_Ir8 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\BR_IP_DNoise_Irr8_perms1000_CI95.mat');
    % Distance Noise, 10mW Source, 5mW Threshold
    FRC_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\FR_IMA_DNoise_Irr10_perms1000_CI95.mat');
    FRP_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\FR_IP_DNoise_Irr10_perms1000_CI95.mat');
    DRC_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\DR_IMA_DNoise_Irr10_perms1000_CI95.mat');
    DRP_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\DR_IP_DNoise_Irr10_perms1000_CI95.mat');
    BRC_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\BR_IMA_DNoise_Irr10_perms1000_CI95.mat');
    BRP_Ir10 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\BR_IP_DNoise_Irr10_perms1000_CI95.mat');

    % Plotting Dist Noise Duration Shuffle Comparisons
    pStruct.yUp = [4.5 15];
    pStruct.yDn = [0 5];

    % FR IMA
    pStruct.WFC = 'FRC';
    pStruct.legendCell = {'No Scattering','Source 8mW','Source 10mW','95% CI'};
    pStruct.legendLoc = 'northeast';
    [caa,daa] = plotSumComp(FRC_NDN.FR_IMA_shufs,FRC_Ir8.FR_IMA_shufs,FRC_Ir10.FR_IMA_shufs,pStruct);

    % DR IMA
    pStruct.WFC = 'DRC';
    pStruct.legendCell = 0;
    [cab,dab] = plotSumComp(DRC_NDN.DR_IMA_shufs,DRC_Ir8.DR_IMA_shufs,DRC_Ir10.DR_IMA_shufs,pStruct);
    
    % BR IMA
    pStruct.WFC = 'BRC';
    [cac,dac] = plotSumComp(BRC_NDN.BR_IMA_shufs,BRC_Ir8.BR_IMA_shufs,BRC_Ir10.BR_IMA_shufs,pStruct);

    % FR IP
    pStruct.WFC = 'FRP';
    [cad,dad] = plotSumComp(FRP_NDN.FR_IP_shufs,FRP_Ir8.FR_IP_shufs,FRP_Ir10.FR_IP_shufs,pStruct);

    % DR IP
    pStruct.WFC = 'DRP';
    [cae,dae] = plotSumComp(DRP_NDN.DR_IP_shufs,DRP_Ir8.DR_IP_shufs,DRP_Ir10.DR_IP_shufs,pStruct);
    
    % BR IP
    pStruct.WFC = 'BRP';
    [caf,daf] = plotSumComp(BRP_NDN.BR_IP_shufs,BRP_Ir8.BR_IP_shufs,BRP_Ir10.BR_IP_shufs,pStruct);

    % Save
    saveDir = 'F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\DistNoise_searches\';
    saveas(caa,[saveDir,'FR_IMA_shufflesComp'],'fig')
    saveas(caa,[saveDir,'FR_IMA_shufflesComp'],'svg')
    saveas(cab,[saveDir,'DR_IMA_shufflesComp'],'fig')
    saveas(cab,[saveDir,'DR_IMA_shufflesComp'],'svg')
    saveas(cac,[saveDir,'BR_IMA_shufflesComp'],'fig')
    saveas(cac,[saveDir,'BR_IMA_shufflesComp'],'svg')
    saveas(cad,[saveDir,'FR_IP_shufflesComp'],'fig')
    saveas(cad,[saveDir,'FR_IP_shufflesComp'],'svg')
    saveas(cae,[saveDir,'DR_IP_shufflesComp'],'fig')
    saveas(cae,[saveDir,'DR_IP_shufflesComp'],'svg')
    saveas(caf,[saveDir,'BR_IP_shufflesComp'],'fig')
    saveas(caf,[saveDir,'BR_IP_shufflesComp'],'svg')
    
    saveas(daa,[saveDir,'FR_IMA_seqShufflesComp'],'fig')
    saveas(daa,[saveDir,'FR_IMA_seqShufflesComp'],'svg')
    saveas(dab,[saveDir,'DR_IMA_seqShufflesComp'],'fig')
    saveas(dab,[saveDir,'DR_IMA_seqShufflesComp'],'svg')
    saveas(dac,[saveDir,'BR_IMA_seqShufflesComp'],'fig')
    saveas(dac,[saveDir,'BR_IMA_seqShufflesComp'],'svg')
    saveas(dad,[saveDir,'FR_IP_seqShufflesComp'],'fig')
    saveas(dad,[saveDir,'FR_IP_seqShufflesComp'],'svg')
    saveas(dae,[saveDir,'DR_IP_seqShufflesComp'],'fig')
    saveas(dae,[saveDir,'DR_IP_seqShufflesComp'],'svg')
    saveas(daf,[saveDir,'BR_IP_seqShufflesComp'],'fig')
    saveas(daf,[saveDir,'BR_IP_seqShufflesComp'],'svg')
    
    close all
    
end


%% Opsin Expression Noise
if opsinNoiseFlag == 1

    % Load Variables
    % No Opsin Noise
    FRC_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IMA_1000_CI95.mat');
    FRP_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IP_1000_CI95.mat');
    DRC_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IMA_1000_CI95.mat');
    DRP_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IP_1000_CI95.mat');
    BRC_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IMA_1000_CI95.mat');
    BRP_NON = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IP_1000_CI95.mat');
    % Opsin Noise, 0.05
    FRC_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\FR_IMA_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    FRP_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\FR_IP_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    DRC_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\DR_IMA_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    DRP_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\DR_IP_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    BRC_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\BR_IMA_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    BRP_Sg05 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\BR_IP_ChR2Noise_sigma0_05_perms1000_CI95.mat');
    % Distance Noise, 10mW Source, 5mW Threshold
    FRC_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\FR_IMA_ChR2Noise_sigma0_1_perms1000_CI95.mat');
    FRP_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\FR_IP_ChR2Noise_sigma0_1_perms1000_CI95.mat');
    DRC_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\DR_IMA_ChR2Noise_sigma0_1_perms1000_CI95.mat');
    DRP_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\DR_IP_ChR2Noise_sigma0_1_perms1000_CI95.mat');
    BRC_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\BR_IMA_ChR2Noise_sigma0_1_perms1000_CI95.mat');
    BRP_Sg1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\BR_IP_ChR2Noise_sigma0_1_perms1000_CI95.mat');

    % Plotting Dist Noise Duration Shuffle Comparisons
    pStruct.yUp = [4 15];
    pStruct.yDn = [0 5];

    % FR IMA
    pStruct.WFC = 'FRC';
    pStruct.legendCell = {'No Opsin Noise','Noise \sigma 0.05','Noise \sigma 0.10','95% CI'};
    pStruct.legendLoc = 'northeast';
    [eaa,faa] = plotSumComp(FRC_NON.FR_IMA_shufs,FRC_Sg05.FR_IMA_shufs,FRC_Sg1.FR_IMA_shufs,pStruct);

    % DR IMA
    pStruct.WFC = 'DRC';
    pStruct.legendCell = 0;
    [eab,fab] = plotSumComp(DRC_NON.DR_IMA_shufs,DRC_Sg05.DR_IMA_shufs,DRC_Sg1.DR_IMA_shufs,pStruct);
    
    % BR IMA
    pStruct.WFC = 'BRC';
    [eac,fac] = plotSumComp(BRC_NON.BR_IMA_shufs,BRC_Sg05.BR_IMA_shufs,BRC_Sg1.BR_IMA_shufs,pStruct);

    % FR IP
    pStruct.WFC = 'FRP';
    [ead,fad] = plotSumComp(FRP_NON.FR_IP_shufs,FRP_Sg05.FR_IP_shufs,FRP_Sg1.FR_IP_shufs,pStruct);

    % DR IP
    pStruct.WFC = 'DRP';
    [eae,fae] = plotSumComp(DRP_NON.DR_IP_shufs,DRP_Sg05.DR_IP_shufs,DRP_Sg1.DR_IP_shufs,pStruct);
    
    % BR IP
    pStruct.WFC = 'BRP';
    [eaf,faf] = plotSumComp(BRP_NON.BR_IP_shufs,BRP_Sg05.BR_IP_shufs,BRP_Sg1.BR_IP_shufs,pStruct);

    % Save
    saveDir = 'F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\ChR2Noise_searches\';
    saveas(eaa,[saveDir,'FR_IMA_shufflesComp'],'fig')
    saveas(eaa,[saveDir,'FR_IMA_shufflesComp'],'svg')
    saveas(eab,[saveDir,'DR_IMA_shufflesComp'],'fig')
    saveas(eab,[saveDir,'DR_IMA_shufflesComp'],'svg')
    saveas(eac,[saveDir,'BR_IMA_shufflesComp'],'fig')
    saveas(eac,[saveDir,'BR_IMA_shufflesComp'],'svg')
    saveas(ead,[saveDir,'FR_IP_shufflesComp'],'fig')
    saveas(ead,[saveDir,'FR_IP_shufflesComp'],'svg')
    saveas(eae,[saveDir,'DR_IP_shufflesComp'],'fig')
    saveas(eae,[saveDir,'DR_IP_shufflesComp'],'svg')
    saveas(eaf,[saveDir,'BR_IP_shufflesComp'],'fig')
    saveas(eaf,[saveDir,'BR_IP_shufflesComp'],'svg')
    
    saveas(faa,[saveDir,'FR_IMA_seqShufflesComp'],'fig')
    saveas(faa,[saveDir,'FR_IMA_seqShufflesComp'],'svg')
    saveas(fab,[saveDir,'DR_IMA_seqShufflesComp'],'fig')
    saveas(fab,[saveDir,'DR_IMA_seqShufflesComp'],'svg')
    saveas(fac,[saveDir,'BR_IMA_seqShufflesComp'],'fig')
    saveas(fac,[saveDir,'BR_IMA_seqShufflesComp'],'svg')
    saveas(fad,[saveDir,'FR_IP_seqShufflesComp'],'fig')
    saveas(fad,[saveDir,'FR_IP_seqShufflesComp'],'svg')
    saveas(fae,[saveDir,'DR_IP_seqShufflesComp'],'fig')
    saveas(fae,[saveDir,'DR_IP_seqShufflesComp'],'svg')
    saveas(faf,[saveDir,'BR_IP_seqShufflesComp'],'fig')
    saveas(faf,[saveDir,'BR_IP_seqShufflesComp'],'svg')
    
    close all

end

%% Voltage Noise
if VNoiseFlag == 1

    % Load Variables
    % No Distance Noise
    FRC_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IMA_1000_CI95.mat');
    FRP_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\FR_IP_1000_CI95.mat');
    DRC_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IMA_1000_CI95.mat');
    DRP_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\DR_IP_1000_CI95.mat');
    BRC_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IMA_1000_CI95.mat');
    BRP_NVN = load('F:\Research\Manuscripts\Computational Waveform\Final Code_V1\ComputationalWaveformMemory\Replay Extension Final\Variables\BR_IP_1000_CI95.mat');
    % Distance Noise, 8mW Source, 5mW Threshold
    FRC_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\FR_IMA_NAmp0_2_perms1000_CI95.mat');
    FRP_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\FR_IP_NAmp0_2_perms1000_CI95.mat');
    DRC_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\DR_IMA_NAmp0_2_perms1000_CI95.mat');
    DRP_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\DR_IP_NAmp0_2_perms1000_CI95.mat');
    BRC_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\BR_IMA_NAmp0_2_perms1000_CI95.mat');
    BRP_Amp02 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\BR_IP_NAmp0_2_perms1000_CI95.mat');
    % Distance Noise, 10mW Source, 5mW Threshold
    FRC_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\FR_IMA_NAmp1_perms1000_CI95.mat');
    FRP_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\FR_IP_NAmp1_perms1000_CI95.mat');
    DRC_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\DR_IMA_NAmp1_perms1000_CI95.mat');
    DRP_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\DR_IP_NAmp1_perms1000_CI95.mat');
    BRC_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\BR_IMA_NAmp1_perms1000_CI95.mat');
    BRP_Amp1 = load('F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\BR_IP_NAmp1_perms1000_CI95.mat');

    % Plotting Dist Noise Duration Shuffle Comparisons
    pStruct.yUp = [3 15];
    pStruct.yDn = [0 5];

    % FR IMA
    pStruct.WFC = 'FRC';
    pStruct.legendCell = {'No Voltage Noise', 'V Noise Amp 0.20','V Noise Amp 1.00','95% CI'};
    pStruct.legendLoc = 'northeast';
    [gaa,haa] = plotSumComp(FRC_NVN.FR_IMA_shufs,FRC_Amp02.FR_IMA_shufs,FRC_Amp1.FR_IMA_shufs,pStruct);

    % DR IMA
    pStruct.WFC = 'DRC';
    pStruct.legendCell = 0;
    [gab,hab] = plotSumComp(DRC_NVN.DR_IMA_shufs,DRC_Amp02.DR_IMA_shufs,DRC_Amp1.DR_IMA_shufs,pStruct);
    
    % BR IMA
    pStruct.WFC = 'BRC';
    [gac,hac] = plotSumComp(BRC_NVN.BR_IMA_shufs,BRC_Amp02.BR_IMA_shufs,BRC_Amp1.BR_IMA_shufs,pStruct);

    % FR IP
    pStruct.WFC = 'FRP';
    [gad,had] = plotSumComp(FRP_NVN.FR_IP_shufs,FRP_Amp02.FR_IP_shufs,FRP_Amp1.FR_IP_shufs,pStruct);

    % DR IP
    pStruct.WFC = 'DRP';
    [gae,hae] = plotSumComp(DRP_NVN.DR_IP_shufs,DRP_Amp02.DR_IP_shufs,DRP_Amp1.DR_IP_shufs,pStruct);
    
    % BR IP
    pStruct.WFC = 'BRP';
    [gaf,haf] = plotSumComp(BRP_NVN.BR_IP_shufs,BRP_Amp02.BR_IP_shufs,BRP_Amp1.BR_IP_shufs,pStruct);

    % Save
    saveDir = 'F:\Research\Code\CA3 Region Code\Ripple Extension Model\RipExtend Figures V2\RipExtend_Noise\VoltageNoise_searches\';
    saveas(gaa,[saveDir,'FR_IMA_shufflesComp'],'fig')
    saveas(gaa,[saveDir,'FR_IMA_shufflesComp'],'svg')
    saveas(gab,[saveDir,'DR_IMA_shufflesComp'],'fig')
    saveas(gab,[saveDir,'DR_IMA_shufflesComp'],'svg')
    saveas(gac,[saveDir,'BR_IMA_shufflesComp'],'fig')
    saveas(gac,[saveDir,'BR_IMA_shufflesComp'],'svg')
    saveas(gad,[saveDir,'FR_IP_shufflesComp'],'fig')
    saveas(gad,[saveDir,'FR_IP_shufflesComp'],'svg')
    saveas(gae,[saveDir,'DR_IP_shufflesComp'],'fig')
    saveas(gae,[saveDir,'DR_IP_shufflesComp'],'svg')
    saveas(gaf,[saveDir,'BR_IP_shufflesComp'],'fig')
    saveas(gaf,[saveDir,'BR_IP_shufflesComp'],'svg')
    
    saveas(haa,[saveDir,'FR_IMA_seqShufflesComp'],'fig')
    saveas(haa,[saveDir,'FR_IMA_seqShufflesComp'],'svg')
    saveas(hab,[saveDir,'DR_IMA_seqShufflesComp'],'fig')
    saveas(hab,[saveDir,'DR_IMA_seqShufflesComp'],'svg')
    saveas(hac,[saveDir,'BR_IMA_seqShufflesComp'],'fig')
    saveas(hac,[saveDir,'BR_IMA_seqShufflesComp'],'svg')
    saveas(had,[saveDir,'FR_IP_seqShufflesComp'],'fig')
    saveas(had,[saveDir,'FR_IP_seqShufflesComp'],'svg')
    saveas(hae,[saveDir,'DR_IP_seqShufflesComp'],'fig')
    saveas(hae,[saveDir,'DR_IP_seqShufflesComp'],'svg')
    saveas(haf,[saveDir,'BR_IP_seqShufflesComp'],'fig')
    saveas(haf,[saveDir,'BR_IP_seqShufflesComp'],'svg')
    
    close all

end
