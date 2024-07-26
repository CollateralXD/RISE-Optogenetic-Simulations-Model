%%% plot learn model summary plots
% Plots comparison of shuffles for the ripple extension model
% LKW 9/21/2021
set(0,'DefaultLineLineWidth',2)

adaptFlag = 1; 
distNoiseFlag = 0;
opsinNoiseFlag = 0;
VNoiseFlag = 0;

pStruct.nRamps = 21;
pStruct.nPs = 101;

%% For Adaptation Sum Plots
if adaptFlag == 1
    
    % Load Variables
    % No Adaptation Ww = 0.029
    FRC_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\FR_IMA_learnDur_80_perms500_CI95.mat');
    DRC_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\DR_IMA_learnDur_80_perms500_CI95.mat');
    BRC_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\BR_IMA_learnDur_80_perms500_CI95.mat');
    FRP_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\FR_IP_learnDur_80_perms500_CI95.mat');
    DRP_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\DR_IP_learnDur_80_perms500_CI95.mat');
    BRP_NA_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\ripLearn_NoAdapt_rerun\BR_IP_learnDur_80_perms500_CI95.mat');
    % Adaptation Ww = 0.029
    FRC_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\FR_IMA_WwMax0_031_learnDur80_perms500_CI95.mat');
    DRC_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\DR_IMA_WwMax0_031_learnDur80_perms500_CI95.mat');
    BRC_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\BR_IMA_WwMax0_031_learnDur80_perms500_CI95.mat');
    FRP_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\FR_IP_WwMax0_031_learnDur80_perms500_CI95.mat');
    DRP_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\DR_IP_WwMax0_031_learnDur80_perms500_CI95.mat');
    BRP_A_W031 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\BR_IP_WwMax0_031_learnDur80_perms500_CI95.mat');
    %Adaptation Ww = 0.032
    FRC_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\FR_IMA_WwMax0_035_learnDur80_perms500_CI95.mat');
    DRC_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\DR_IMA_WwMax0_035_learnDur80_perms500_CI95.mat');
    BRC_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\BR_IMA_WwMax0_035_learnDur80_perms500_CI95.mat');
    FRP_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\FR_IP_WwMax0_035_learnDur80_perms500_CI95.mat');
    DRP_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\DR_IP_WwMax0_035_learnDur80_perms500_CI95.mat');
    BRP_A_W035 = load('F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\BR_IP_WwMax0_035_learnDur80_perms500_CI95.mat');
    
    % Plotting Pulse Duration Shuffle Comparisons
    pStruct.yUp = [4 15];
    pStruct.yDn = [0 0];
    
    % FR IMA
    pStruct.WFC = 'FRC';
    pStruct.legendCell = {'No Adapt Ww=0.029','Adapt Ww=0.031','Adapt Ww=0.035','95% CI'};
    pStruct.legendLoc = 'northeast';
    [aaa,baa] = plotSumComp_learnfun(FRC_NA_W031.FR_IMA_shufs,FRC_A_W031.FR_IMA_shufs,FRC_A_W035.FR_IMA_shufs,pStruct);

    % DR IMA
    pStruct.WFC = 'DRC';
    pStruct.legendCell = 0;
    [aab,bab] = plotSumComp_learnfun(DRC_NA_W031.DR_IMA_shufs,DRC_A_W031.DR_IMA_shufs,DRC_A_W035.DR_IMA_shufs,pStruct);
    
    % BR IMA
    pStruct.WFC = 'BRC';
    [aac,bac] = plotSumComp_learnfun(BRC_NA_W031.BR_IMA_shufs,BRC_A_W031.BR_IMA_shufs,BRC_A_W035.BR_IMA_shufs,pStruct);

    % FR IP
    pStruct.WFC = 'FRP';
    [aad,bad] = plotSumComp_learnfun(FRP_NA_W031.FR_IP_shufs,FRP_A_W031.FR_IP_shufs,FRP_A_W035.FR_IP_shufs,pStruct);

    % DR IP
    pStruct.WFC = 'DRP';
    [aae,bae] = plotSumComp_learnfun(DRP_NA_W031.DR_IP_shufs,DRP_A_W031.DR_IP_shufs,DRP_A_W035.DR_IP_shufs,pStruct);
    
    % BR IP
    pStruct.WFC = 'BRP';
    [aaf,baf] = plotSumComp_learnfun(BRP_NA_W031.BR_IP_shufs,BRP_A_W031.BR_IP_shufs,BRP_A_W035.BR_IP_shufs,pStruct);
    
%     % Save
%     saveDir = 'F:\Research\Code\CA3 Region Code\Pattern Learning Model\LearningFigures_v2\learnAdapt_Figs\';
%     saveas(aaa,[saveDir,'FR_IMA_shufflesComp'],'fig')
%     saveas(aaa,[saveDir,'FR_IMA_shufflesComp'],'svg')
%     saveas(aab,[saveDir,'DR_IMA_shufflesComp'],'fig')
%     saveas(aab,[saveDir,'DR_IMA_shufflesComp'],'svg')
%     saveas(aac,[saveDir,'BR_IMA_shufflesComp'],'fig')
%     saveas(aac,[saveDir,'BR_IMA_shufflesComp'],'svg')
%     saveas(aad,[saveDir,'FR_IP_shufflesComp'],'fig')
%     saveas(aad,[saveDir,'FR_IP_shufflesComp'],'svg')
%     saveas(aae,[saveDir,'DR_IP_shufflesComp'],'fig')
%     saveas(aae,[saveDir,'DR_IP_shufflesComp'],'svg')
%     saveas(aaf,[saveDir,'BR_IP_shufflesComp'],'fig')
%     saveas(aaf,[saveDir,'BR_IP_shufflesComp'],'svg')
%     
%     saveas(baa,[saveDir,'FR_IMA_seqShufflesComp'],'fig')
%     saveas(baa,[saveDir,'FR_IMA_seqShufflesComp'],'svg')
%     saveas(bab,[saveDir,'DR_IMA_seqShufflesComp'],'fig')
%     saveas(bab,[saveDir,'DR_IMA_seqShufflesComp'],'svg')
%     saveas(bac,[saveDir,'BR_IMA_seqShufflesComp'],'fig')
%     saveas(bac,[saveDir,'BR_IMA_seqShufflesComp'],'svg')
%     saveas(bad,[saveDir,'FR_IP_seqShufflesComp'],'fig')
%     saveas(bad,[saveDir,'FR_IP_seqShufflesComp'],'svg')
%     saveas(bae,[saveDir,'DR_IP_seqShufflesComp'],'fig')
%     saveas(bae,[saveDir,'DR_IP_seqShufflesComp'],'svg')
%     saveas(baf,[saveDir,'BR_IP_seqShufflesComp'],'fig')
%     saveas(baf,[saveDir,'BR_IP_seqShufflesComp'],'svg')
%     
%     close all
% 
end
