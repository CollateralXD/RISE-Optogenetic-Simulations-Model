%%% Sample Code for Learning Model Search Script
%Relies on ca3Seq_learn_V2.m function, learn_CohenStats.m ca3learn_stats.m and
%the plot_CA3LearnAct.m functions
%Searches ripple prolongation parameter space in dimensions of ramp degree
%and input pulse duration
%Plots various metrics and statistics for evaluating network performance
%under varying waveform shape input
%LKW 8/17/21
%%%

clearvars
close all
load learnControl_adapt_sym_offAxis1_37.mat

pStruct.rampTypeFlag    = 1;    %1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 6 = BR IMA; 7 = BR IP
suppGraphFlag           = 0;    %Plot extra graphs or not
pStruct.simTypeFlag     = 2;    %1 = Linear; 2 = Linear with Adaptation
pStruct.noiseFlag       = 0;    %0 = no noise; 1 = White noise; 2 = 1/f or pink noise
pStruct.wtNoiseFlag     = 0;    %0 = no noise; 1 = Lognormal noise distribution using parameters Mu and Sigma
gutCheckFlag            = 1;    %0 = no gut check; 1 = gut check
saveFlag                = 0;    %0 = no save;      1 = saves
saveDir = 'Your\Dir\Here\';
prmtag = '_altCtlWt1_37_learnDur';            %name of changed parameter

disp(['Ramp Type ', num2str(pStruct.rampTypeFlag),'; Sim Type ', num2str(pStruct.simTypeFlag), ...
    '; Noise type ', num2str(pStruct.noiseFlag), '; Wt Mat Noise type ', num2str(pStruct.wtNoiseFlag)]);

%Setup neuron and weight paramters
pStruct.N               = 15;       %Nodes per region
N                       = 15;       %For gut check run
dt                      = 1;
pStruct.Ww              = 0.035;    %Weight strength pyr to pyr
pStruct.Hh              = 0.05;     %Weight strength IN to Pyr
pStruct.Wh              = 0.05;     %Weight strength pyr to IN
pStruct.HAuto           = 0.003;    %Inhibition self feedback
pStruct.tha             = 4*ones(1,N);
pStruct.thh             = 4*ones(1,N);
pStruct.eta1            = 0.01;     %Decay constant
pStruct.eta2            = 0.01;
pStruct.actThresh       = 10*ones(1,N);
pStruct.e1              = 0.002;     %Presynaptic decay
if pStruct.wtNoiseFlag == 1
    pStruct.e3          = 0.01;     %Learning rate
    pStruct.achLearn    = 0.1;       %Acetylcholine state; 0 = high Ach; 1 = low Ach
else
    pStruct.e3          = 0.001;
    pStruct.achLearn    = 0.1;       %Acetylcholine state; 0 = high Ach; 1 = low Ach
end
pStruct.noiseAmp        = 0.1*pStruct.Ww;          %Amplitude of membrane noise if included.
pStruct.noiseMu         = 1;            %Mean of lognormal wt distro
pStruct.noiseSigma      = 0.4;            %Variance of lognorm wt distro

%Setup input parameters
pStruct.dt              = 1;
pStruct.TLearn          = 1500/dt;     %Time steps, must be even
pStruct.TTest           = 1500/dt;
pStruct.cueN            = 1;
pStruct.Iexcit1         = 0.5;      %Strength of learn pulses (Iso current)
pStruct.Iexcit2         = 1;     	%Strength of test cue pulse
pStruct.testDur         = 20/dt;       %Test phase cue duration
pStruct.learnDur        = 80/dt;       %Learn phase pulse duration
% pStruct.learnOverlap    = 44;       %Learn phase pulse overlap (ms)
pStruct.learnOverlap    = 0.55;      %Learn phase pulse overlap %
pStruct.rampPerc        = 0.25;
pStruct.rampLen1        = round(pStruct.learnDur*pStruct.rampPerc);  %Duration ramp length
% pStruct.rampLen2        = round(pStruct.learnDur*pStruct.rampPerc);
pStruct.onsetDelay      = 50/dt;        %Wait time to ripple start from sim start
pStruct.achTest         = 1;

% Ionic Currents and related parameters 
pStruct.mu              = 0.01;      %Ca-dependent K-current
pStruct.gm              = 0.001;     %Gamma; voltage-dependent Ca-currents
pStruct.om              = 0.001;     %Omega; constant for diffusion of intracellular Ca
pStruct.thc             = 4*ones(1,N);

W                       = zeros(N,N,pStruct.TLearn);           %Init wt mat
AH                      = zeros(N,N);           %Wt mat of a to h (feedforward activation of IN)
H                       = zeros(N,N)*pStruct.Hh/10;           %Weight of h to a (feedback inhibition)
ahGain                  = 5;
hGain                   = 5;
for i = 1:N         %Build weight mats
    H(i,i)  = pStruct.Hh;   %Direct feedback IN to Pyr
    AH(i,i) = pStruct.Wh;   %Direct excitation pyr to IN
    if i <= N - 1   %Forward 1
%         H(i,i+1)  = pStruct.Hh/hGain;   %IN  2 Pyr
%         AH(i,i+1) = pStruct.Wh/ahGain;   %Pyr 2 IN
    end
    if i <= N - 2   %Forward 2
%         H(i,i+2)  = pStruct.Hh/hGain;   %IN  2 Pyr
%         AH(i,i+2) = pStruct.Wh/ahGain;   %Pyr 2 IN
    end
    if i <= N - 3   %Forward 3
%         H(i,i+3)  = pStruct.Hh/4;   %IN  2 Pyr
%         AH(i,i+3) = pStruct.Wh/4;   %Pyr 2 IN
    end
    if i > 1        %Back 1
%         H(i,i-1)  = pStruct.Hh/hGain;   %IN  2 Pyr
%         AH(i,i-1) = pStruct.Wh/ahGain;   %Pyr 2 IN
    end
    if i > 2        %Back 2   
%         H(i,i-2)  = pStruct.Hh/hGain;   %IN  2 Pyr
%         AH(i,i-2) = pStruct.Wh/ahGain;   %Pyr 2 IN
    end
    if i > 3        %Back 3
%         H(i,i-3)  = pStruct.Hh/1;   %IN  2 Pyr
%         AH(i,i-3) = pStruct.Wh/4;   %Pyr 2 IN
    end
end

pStruct.W = W; pStruct.H = H; pStruct.AH = AH;  %Add wtmats to pStruct

%Set up search variables and parameters
pFloor                  = 0;
pLim                    = 1;        %1 for % learnOverlap
nPs                     = 61;
% pLim                    = 150;  %150 for learnDur
% nPs                     = 61;
nRamps                  = 21;
pVect                   = linspace(pFloor,pLim,nPs);    %Vector of parameter values
if mod(pStruct.rampTypeFlag,2) == 1; rampLim = 1; else; rampLim = 0.5; end %Set ramp limit to 1 or 0.5 depending on flag
rampPercs = linspace(0,rampLim,nRamps);

outputs.diffs       = cell(nRamps,nPs);
outputs.INdiffs     = cell(nRamps,nPs);
outputs.testPks     = cell(nRamps,nPs);
outputs.testLocs    = cell(nRamps,nPs);
outputs.Ws          = cell(nRamps,nPs);     %Storage of final Wt. Mat snapshot for each simulation
outputs.Ds          = zeros(nRamps,nPs);

%% Gut check single above parameters
if gutCheckFlag == 1
    [actCell,hactCell,inCell]   = ca3Seq_learn_V2(pStruct);    %Using ramp defined above
    aLearn = actCell{1}; aTest  = actCell{2};
    hLearn = hactCell{1};hTest  = hactCell{2};
    ALearn = inCell{1};  ATest  = inCell{2};
    % dGut = ripExtend_CohenStats(actCell,hactCell,pStruct);
    
    % Calculate Peaks
    pksLearn = zeros(N,1); locsLearn = zeros(N,1);
    pksTest = zeros(N,1); locsTest = zeros(N,1);
    
    for  i = 1:N
        [pksTmp,locsTmp] = findpeaks(aLearn(i,:));
        if ~isempty(locsTmp); pksLearn(i) = pksTmp(1); locsLearn(i) = locsTmp(1); end
        [pksTmp,locsTmp] = findpeaks(aTest(i,:));
        if ~isempty(locsTmp); pksTest(i) = pksTmp(1); locsTest(i) = locsTmp(1); end
    end
    pksTest = pksTest(pksTest>0); locsTest = locsTest(locsTest>0);
    locsCell = {locsLearn,locsTest};
    pksCell = {pksLearn,pksTest};
    
    %Threshold method
    actThresh = 10*ones(1,N); % tha;
    tttCell = {};
    
    for i = 1:2
        pksTmp = pksCell{i};        %All the peaks for that sim (e.g. ramp)
        actTmp = actCell{i};        %Pyramidal activation for that sim e.g. NxT
        ttt = [];
        for j = 1:numel(pksTmp)     %For each peak in the sim
            indtmp = find(actTmp(j,:)>actThresh(j),1);
            if ~isempty(indtmp)
                ttt = [ttt,indtmp];
            end
        end
        tttCell(i) = {ttt};
    end
    
    [thaPlot,thax1,thax2] = plotCA3Ripples_learn(actCell,hactCell,inCell,locsCell,pksCell,pStruct);
    set(thaPlot,'defaultLegendAutoUpdate','off');
    axCell = {thax1,thax2};
    thaLength = 20;
    for i = 1:2
        axes(axCell{i});
        tmpTTT = tttCell{i};
        plot([0,999],[10,10],'k--','LineWidth',1)
        for j = 1:length(tmpTTT)
            %         plot([tmpTTT(j)-thaLength,tmpTTT(j)+thaLength],[tha(j),tha(j)],'k--','LineWidth',1)
            scatter(tmpTTT(j),actThresh(j),'k^','filled')
        end
    end
%     [rampINplot,rampINax1,rampINax2] = plotCA3INs(actCell,hactCell,inCell,pStruct.cueN,1);
%     [controlINplot,controlINax1,controlINax2] = plotCA3INs(actCell,hactCell,inCell,pStruct.cueN,3);
end

%% Run Block Parameter Search

for j = 1:nPs
    pStruct.learnOverlap = pVect(j);    %2nd search parameter: learn OL
%     pStruct.learnDur = pVect(j);        %2nd search parameter: learn pulse dur
%     pStruct.noiseAmp = pVect(j)*pStruct.Ww;
    for i = 1:nRamps
        pStruct.tmpRamp = rampPercs(i);                             %Define ramp percentage
        pStruct.rampLen1= round(pStruct.tmpRamp*pStruct.learnDur);   %Define ramp length
                
        %Run function and calculate stats
        [actTmp,hactTmp,inTmp,outputs.Ws{i,j}] = ca3Seq_learn_V2(pStruct);
        if max(actTmp{2},[],2) > 100
            outputs.testN(i,j) = 0;
            outputs.Ds(i,j) = NaN;
        else
            [outputs.testPks{i,j},outputs.testLocs{i,j},outputs.diffs{i,j},outputs.INDiffs{i,j}] = ca3learn_stats(pStruct,actTmp,hactTmp);
            outputs.testN(i,j) = numel(outputs.testPks{i,j});
            [outputs.Ds(i,j),~,~] = learn_CohenStats({actTmp{2},aControl},NaN,pStruct);
        end
        %Calculate other outputs numbers
        if ~isempty(outputs.diffs{i,j})
            outputs.meanIThI(i,j) = mean(outputs.diffs{i,j});
%             outputs.meanINIThI(i,j) = mean(outputs.INDiffs{i,j});
%             outputs.sdIThI(i,j) = std(outputs.diffs{i,j});
%             outputs.sdINIThI(i,j) = std(outputs.INDiffs{i,j});
        else
            outputs.meanIThI(i,j) = 0;
%             outputs.meanINIThI(i,j) = 0;
        end
    end
end

%% Clean outputs and find lowest OL duration of full replay

% Clean ds output
outputs.nanMaxDs = outputs.Ds;
% outputs.cleanDs = outputs.Ds .* nan((outputs.Ds >= 0));
if pStruct.wtNoiseFlag == 1   %Use 99th percentile instead of data max
    linearDs = abs(reshape(outputs.nanMaxDs,[1 numel(outputs.nanMaxDs)]));
    tmpMax = prctile(linearDs,90)
else                        %Use maximum, non Inf simulation value
    tmpMax = max(max(abs(~isinf(outputs.nanMaxDs).*outputs.nanMaxDs)));
end
% tmpMax = max(max(outputs.nanMaxDs));
outputs.nanMaxDs(isnan(outputs.nanMaxDs)) = tmpMax;  %Convert any nans to max effect size for plotting purposes
outputs.nanZedDs = outputs.Ds;
outputs.nanZedDs(isnan(outputs.nanZedDs)) = 0;      %Convert any nans to 0s for shuffle purposes

bestSeqInds = bestSeqN(flipud(outputs.testN'),pVect,nRamps);
minDs   = min(flipud(abs(outputs.nanMaxDs')));

cbMax = tmpMax;
% cbMax = 4;
%% Plotting

%HeatMap of Pyr temporal disruption (Effect Size)
figure; % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]);
aaa = gca;
imagesc(flipud(abs(outputs.nanMaxDs')),[0,cbMax]);
aaa.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaa.XTick = 1:round(nRamps-1)/5:nRamps;
if pLim == 1; yticklabels(linspace(pLim*100,0,6)); else yticklabels(linspace(pLim,0,6)); end
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('OL Duration (ms)'); 
colormap hot; axis square; hcb = colorbar;
% hcb.Label.String = "Cohen's d";
set(aaa,'FontSize',24,'fontname','times')

%Plot map of sequence length
seqPlot = figure;  aaz = gca;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
imagesc(flipud(outputs.testN'),[0,N]);
aaz.YTick = 1:round((nPs-1)/5):nPs; 
aaz.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
if pLim == 1; yticklabels(linspace(pLim*100,0,6)); else yticklabels(linspace(pLim,0,6)); end
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Learn Overlap (ms)'); 
colormap parula; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Sequence Replay Length; Forward Ramp Learn");
% else; title("Sequence Replay Length; Double Ramp Learn"); end
% hcb.Label.String = "# Suprathreshold Test Nodes";
set(aaz,'FontSize',24,'fontname','times')

%% Shuffles
permN = 500;
z = 1.96;   %confidence level 95% = 1.96, 99% = 2.58
pyrIThIShuf = bootCols(flipud(outputs.meanIThI'),permN);        %Pyr IThI shuffle distribution
actNShuf = bootCols(flipud(outputs.testN'),permN);              %Pyr Seq length shuffle
rampDShuf = bootColDs(flipud(abs(outputs.nanZedDs')),permN);    %2ndary parameter shuffle for columns of ramp %
prmDShuf = bootColDs(flipud(abs(outputs.nanZedDs)),permN);      %Ramp % shuffle for columsn of 2ndary parameter
prmNShuf = bootCols(flipud(outputs.testN),permN);
shufs.meanPyrIThI = mean(pyrIThIShuf);              %Shuffle means
shufs.meanActN = mean(actNShuf);
shufs.meanRampD = mean(rampDShuf);
shufs.meanPrmD = mean(prmDShuf);
shufs.meanPrmN = mean(prmNShuf);
shufs.sdPyrIThI = std(pyrIThIShuf);                 %Shuffle SDs
shufs.sdActN = std(actNShuf);
shufs.sdRampD = std(rampDShuf);
shufs.sdPrmD = std(prmDShuf);
shufs.sdPrmN = std(prmNShuf);

% CIs
shufs.upCIIThI = shufs.meanPyrIThI + z*shufs.sdPyrIThI/sqrt(permN);
shufs.dnCIIThI = shufs.meanPyrIThI - z*shufs.sdPyrIThI/sqrt(permN);
shufs.upCIN = shufs.meanActN + z*shufs.sdActN/sqrt(permN);
shufs.dnCIN = shufs.meanActN - z*shufs.sdActN/sqrt(permN);
shufs.upCId = shufs.meanRampD + z*shufs.sdRampD/sqrt(permN);
shufs.dnCId = shufs.meanRampD - z*shufs.sdRampD/sqrt(permN);
shufs.upCIPrmD = shufs.meanPrmD + z*shufs.sdPrmD/sqrt(permN);
shufs.dnCIPrmD = shufs.meanPrmD - z*shufs.sdPrmD/sqrt(permN);
shufs.upCIPrmN = shufs.meanPrmN + z*shufs.sdPrmN/sqrt(permN);
shufs.dnCIPrmN = shufs.meanPrmN - z*shufs.sdPrmN/sqrt(permN);

%% Supplementary Graphs
% % Shuffles subplots
% % figure(); % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); 
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% subplot(2,2,1); aab = gca;
% tmpMax = max(max(abs(outputs.meanIThI)));
% axis square; hold on;
% plot(rampPercs*100,shufs.meanPyrIThI,'k',rampPercs*100,upCIIThI,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCIIThI,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean IThI (ms)")
% legend('OL Shuffle','99% CI');
% ylim([0 tmpMax]); aab.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,2); aab = gca;
% tmpMax = max(max(abs(outputs.meanIThI)));
% axis square; hold on;
% plot(rampPercs*100,meanDShuf,'k',rampPercs*100,upCId,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCId,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean Effect Size (ms)")
% ylim([0 3]);
% aab.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100);
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,3); aac = gca;
% tmpMax = max(max(abs(outputs.testN)));
% axis square; hold on;
% plot(rampPercs*100,meanActNShuf,'k',rampPercs*100,upCIN,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCIN,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean Suprathreshold Nodes")
% legend('OL Shuffle','99% CI');
% ylim([0 tmpMax]); aac.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,4); aad = gca;
% axis square; hold on;
% plot(rampPercs*100, bestSeqInds,'k-')
% xlabel('Ramp Percentage'); ylabel('OL dur of Best Seq Length');
% ylim([-10 80]);
% set(gca,'FontSize',24,'fontname','times')

% % Plot map of pyramidal mean IthI
% pyrPlot = figure;  aaa = gca;
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% imagesc(flipud(outputs.meanIThI'),[0,max(max(outputs.meanIThI))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs; 
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% % yticklabels(linspace(pLim,pFloor,6));
% % xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% yticklabels(linspace(100,0,6))
% xticklabels(0:20:100);
% % xlabel('Ramp Percentage'); ylabel('Learn Overlap (ms)'); 
% colormap jet; axis square; hcb = colorbar;
% % if mod(pStruct.rampTypeFlag,2) == 1; title("Replay Mean Pyr IThI; Forward Ramp Learn");
% % else; title("Replay Mean Pyr IThI; Double Ramp Learn"); end
% % hcb.Label.String = "Mean Pyr IThI";
% set(aaa,'FontSize',24,'fontname','times')

% % Plot map of IN mean IthI
% INPlot = figure; aaa = gca;
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% imagesc(flipud(outputs.meanINIThI'),[0,max(max(outputs.meanINIThI))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs; 
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% % xlabel('Ramp Percentage'); 
% % ylabel('Learn Overlap (ms)'); 
% colormap jet; axis square; hcb = colorbar;
% % if mod(pStruct.rampTypeFlag,2) == 1; title("Replay Mean IN IThI; Forward Ramp Learn");
% % else; title("Replay Mean IN IThI; Double Ramp Learn"); end
% % hcb.Label.String = "Mean IN IThI";
% set(aaa,'FontSize',24,'fontname','times')

% %Plot Weights
% figure(); set(gcf, 'Colormap', jet) % tightfig();
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.5, 0.75, 0.5]);
% subplot(1,3,1)
% imagesc(W); title('Pyr-Pyr Weight (W)');    axis square
% xlabel('Post-Synaptic Pyr'); ylabel('Pre-Synaptic Pyr')
% set(gca,'FontSize',20)
% subplot(1,3,2)
% imagesc(AH); title("Pyr-IN Weight (W')");   axis square
% xlabel('Post-Synaptic IN'); ylabel('Pre-Synaptic Pyr')
% set(gca,'FontSize',20)
% subplot(1,3,3)
% imagesc(H); title('IN-Pyr Weight (H)');     axis square
% xlabel('Post-Synaptic Pyr'); ylabel('Pre-Synaptic IN')
% set(gca,'FontSize',20)

% % Afferent 3D plot code
% sqGhost = zeros(1,pStruct.TLearn);
% sqGhost(pStruct.onsetDelay+1:pStruct.onsetDelay+pStruct.learnDur) = pStruct.Iexcit1;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.35, 0.45]); gaa = gca;
% w3 = waterfall(1:pStruct.TLearn,1:N,ALearn); hold on;
% w3.LineWidth = 3;
% w3Sq = waterfall(1:pStruct.TLearn,1,sqGhost);
% w3Sq.LineStyle = '--'; w3Sq.EdgeColor = 'r'; w3Sq.LineWidth = 2;
% gaa.YTick = 1:7:15;
% zlim([0 1]); xlim([0 1200])
% xlabel('Time'); ylabel('Node','position',[-148.7513525348877,6.348030863114843,-0.099984467105344]); zlabel('Amplitude')
% view([-13.2608328513592 45.9028169014085]);
% set(gca,'FontSize',24,'fontname','times')

if suppGraphFlag == 1
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\DRC_sumStats_wtBiasOff_offAxis1_7.mat')
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\DRP_sumStats_wtBiasOff_offAxis1_7.mat')
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\FRC_sumStats_wtBiasOff_offAxis1_7.mat')
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\FRP_sumStats_wtBiasOff_offAxis1_7.mat')
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\BRC_sumStats_wtBiasOff_offAxis1_7.mat')
%     load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\BRP_sumStats_wtBiasOff_offAxis1_7.mat')
    
    % Combined mean IThI x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); saa = gca;
    axis square; hold on;
    plot(FR_IMA_shufs.meanPyrIThI,'r-');
    plot(DR_IMA_shufs.meanPyrIThI,'b-');
    plot(BR_IMA_shufs.meanPyrIThI,'c-');
    plot(FR_IP_shufs.meanPyrIThI,'r--');
    plot(DR_IP_shufs.meanPyrIThI,'b--');
    plot(BR_IP_shufs.meanPyrIThI,'c--');
    patch(patchXs,[FR_IMA_shufs.dnCIIThI,fliplr(FR_IMA_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DR_IMA_shufs.dnCIIThI,fliplr(DR_IMA_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IMA_shufs.dnCIIThI,fliplr(BR_IMA_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FR_IP_shufs.dnCIIThI,fliplr(FR_IP_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DR_IP_shufs.dnCIIThI,fliplr(DR_IP_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IP_shufs.dnCIIThI,fliplr(BR_IP_shufs.upCIIThI)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean IThI')
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([25 35]);xlim([1 nRamps]); saa.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined min D x Ramp shuffle - must load in previous data
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sab = gca;
    axis square; hold on;
    plot(FR_IMA_minDs,'r-');
    plot(DR_IMA_minDs,'b-');
    plot(BR_IMA_minDs,'c-');
    plot(FR_IP_minDs,'r--');
    plot(DR_IP_minDs,'b--');
    plot(BR_IP_minDs,'c--');
    xlabel('Ramp Percentage'); ylabel('Min Effect Size')
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
    ylim([0 1]); xlim([1 nRamps]); sab.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined mean D x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sac = gca;
    axis square; hold on;
    plot(FR_IMA_shufs.meanRampD,'r-');
    plot(DR_IMA_shufs.meanRampD,'b-');
    plot(BR_IMA_shufs.meanRampD,'c-');
    plot(FR_IP_shufs.meanRampD,'r--');
    plot(DR_IP_shufs.meanRampD,'b--');
    plot(BR_IP_shufs.meanRampD,'c--');
    patch(patchXs,[FR_IMA_shufs.dnCId,fliplr(FR_IMA_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DR_IMA_shufs.dnCId,fliplr(DR_IMA_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IMA_shufs.dnCId,fliplr(BR_IMA_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FR_IP_shufs.dnCId,fliplr(FR_IP_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DR_IP_shufs.dnCId,fliplr(DR_IP_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IP_shufs.dnCId,fliplr(BR_IP_shufs.upCId)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 3]); xlim([1 nRamps]); sac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined mean Seq Length x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sad = gca;
    axis square; hold on;
    plot(FR_IMA_shufs.meanActN,'r-');
    plot(DR_IMA_shufs.meanActN,'b-');
    plot(BR_IMA_shufs.meanActN,'c-');
    plot(FR_IP_shufs.meanActN,'r--');
    plot(DR_IP_shufs.meanActN,'b--');
    plot(BR_IP_shufs.meanActN,'c--');
    patch(patchXs,[FR_IMA_shufs.dnCIN,fliplr(FR_IMA_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DR_IMA_shufs.dnCIN,fliplr(DR_IMA_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IMA_shufs.dnCIN,fliplr(BR_IMA_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FR_IP_shufs.dnCIN,fliplr(FR_IP_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DR_IP_shufs.dnCIN,fliplr(DR_IP_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IP_shufs.dnCIN,fliplr(BR_IP_shufs.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean Sequence Length')
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 15]); xlim([1 nRamps]); sad.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')

    % Combined mean D x Parameter shuffle - must load in previous data
    patchXs = [1:nPs,linspace(nPs,1,nPs)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sae = gca;
    axis square; hold on;
    plot(FR_IMA_shufs.meanPrmD,'r-');
    plot(DR_IMA_shufs.meanPrmD,'b-');
    plot(BR_IMA_shufs.meanPrmD,'c-');
    plot(FR_IP_shufs.meanPrmD,'r--');
    plot(DR_IP_shufs.meanPrmD,'b--');
    plot(BR_IP_shufs.meanPrmD,'c--');
    patch(patchXs,[FR_IMA_shufs.dnCIPrmD,fliplr(FR_IMA_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DR_IMA_shufs.dnCIPrmD,fliplr(DR_IMA_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IMA_shufs.dnCIPrmD,fliplr(BR_IMA_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FR_IP_shufs.dnCIPrmD,fliplr(FR_IP_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DR_IP_shufs.dnCIPrmD,fliplr(DR_IP_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IP_shufs.dnCIPrmD,fliplr(BR_IP_shufs.upCIPrmD)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Learn Pulse Length (ms)'); ylabel('Mean Effect Size')
%     legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 5]); xlim([1 nPs]); sae.XTick = linspace(0,nPs,6); xticklabels(0:30:150)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined mean Seq Length x Parameter shuffle - must load in previous data
    patchXs = [1:nPs,linspace(nPs,1,nPs)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); saf = gca;
    axis square; hold on;
    plot(FR_IMA_shufs.meanPrmN,'r-');
    plot(DR_IMA_shufs.meanPrmN,'b-');
    plot(BR_IMA_shufs.meanPrmN,'c-');
    plot(FR_IP_shufs.meanPrmN,'r--');
    plot(DR_IP_shufs.meanPrmN,'b--');
    plot(BR_IP_shufs.meanPrmN,'c--');
    patch(patchXs,[FR_IMA_shufs.dnCIPrmN,fliplr(FR_IMA_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DR_IMA_shufs.dnCIPrmN,fliplr(DR_IMA_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IMA_shufs.dnCIPrmN,fliplr(BR_IMA_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FR_IP_shufs.dnCIPrmN,fliplr(FR_IP_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DR_IP_shufs.dnCIPrmN,fliplr(DR_IP_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BR_IP_shufs.dnCIPrmN,fliplr(BR_IP_shufs.upCIPrmN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Learn Pulse Length (ms)'); ylabel('Mean Sequence Length')
    % legend('FR Equal Max','DR Equal Max','FR Equal Power','DR Equal Power','95% CI','FontSize',16,'location','northwest');
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 15]); xlim([1 nPs]); saf.XTick = linspace(0,nPs,6); xticklabels(0:30:150)
    set(gca,'FontSize',24,'fontname','times')
end

%% Saves
if saveFlag == 1
disp('saving vars and figs')
if z == 1.96
    ciStr = ['_perms',num2str(permN),'_CI95'];
elseif z == 2.58
    ciStr = ['_perms',num2str(permN),'_CI99'];
end
prmtmp = round(pStruct.learnDur);  %value of changed parameter
% prmtmp = round(pStruct.learnDur);
figtag = [prmtag,num2str(prmtmp)]; %name-value string for figure saves

if pStruct.rampTypeFlag == 1
    sBase = ['FR_IMA',figtag,ciStr];
    FR_IMA_minInds  = bestSeqInds;
    FR_IMA_minDs    = minDs;
    FR_IMA_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName,'FR_IMA_minInds','FR_IMA_minDs','FR_IMA_shufs')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 2
    sBase = ['DR_IMA',figtag,ciStr];
    DR_IMA_minInds  = bestSeqInds;
    DR_IMA_minDs    = minDs;
    DR_IMA_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName, 'DR_IMA_minInds','DR_IMA_minDs','DR_IMA_shufs')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 6
    sBase = ['BR_IMA',figtag,ciStr];
    BR_IMA_minInds  = bestSeqInds;
    BR_IMA_minDs    = minDs;
    BR_IMA_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName, 'BR_IMA_minInds','BR_IMA_minDs','BR_IMA_shufs')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 3
    sBase = ['FR_IP',figtag,ciStr];
    FR_IP_minInds  = bestSeqInds;
    FR_IP_minDs    = minDs;
    FR_IP_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName, 'FR_IP_minInds','FR_IP_minDs','FR_IP_shufs')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 4
    sBase = ['DR_IP',figtag,ciStr];
    DR_IP_minInds  = bestSeqInds;
    DR_IP_minDs    = minDs;
    DR_IP_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName, 'DR_IP_minInds','DR_IP_minDs','DR_IP_shufs')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 7
    sBase = ['BR_IP',figtag,ciStr];
    BR_IP_minInds  = bestSeqInds;
    BR_IP_minDs    = minDs;
    BR_IP_shufs    = shufs;
    fName = [saveDir, sBase];
    save(fName, 'BR_IP_minInds','BR_IP_minDs','BR_IP_shufs')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'svg')
end

    if suppGraphFlag == 1
        
        saveas(saa,[saveDir,'allShufs_rampsXmeanIThIs',figtag],'fig')
        saveas(saa,[saveDir,'allShufs_rampsXmeanIThIs',figtag],'png')
        saveas(saa,[saveDir,'allShufs_rampsXmeanIThIs',figtag],'svg')
        saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'fig')
        saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'png')
        saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'svg')
        saveas(sac,[saveDir,'allShufs_rampsXmeanDs',figtag],'fig')
        saveas(sac,[saveDir,'allShufs_rampsXmeanDs',figtag],'png')
        saveas(sac,[saveDir,'allShufs_rampsXmeanDs',figtag],'svg')
        saveas(sad,[saveDir,'allShufs_rampsXmeanNs',figtag],'fig')
        saveas(sad,[saveDir,'allShufs_rampsXmeanNs',figtag],'png')
        saveas(sad,[saveDir,'allShufs_rampsXmeanNs',figtag],'svg')
        saveas(sae,[saveDir,'allShufs_prmsXmeanDs',figtag],'fig')
        saveas(sae,[saveDir,'allShufs_prmsXmeanDs',figtag],'png')
        saveas(sae,[saveDir,'allShufs_prmsXmeanDs',figtag],'svg')
        saveas(saf,[saveDir,'allShufs_prmsXmeanNs',figtag],'fig')
        saveas(saf,[saveDir,'allShufs_prmsXmeanNs',figtag],'png')
        saveas(saf,[saveDir,'allShufs_prmsXmeanNs',figtag],'svg')        
    end
    
disp('Done!')
end

%% Functions

function [shufDistro] = bootCols(sampMat,iters)
%Equal weight for any nonzero, non 1 value, 0 weight for 0s,1s
shufDistro = zeros(iters,size(sampMat,2));
for i = 1:size(sampMat,2)   %For number of columns
    if sum(sampMat(:,i)) <= size(sampMat,1)  %If column is all 1's
        shufDistro(:,i) = mean(sampMat(:,i));
    else
    tmpWts = sampMat(:,i) > 1;
    tmpWts = ones(size(sampMat,1),1).*tmpWts;   %Set wt 0 for any 1 or 0
    shufDistro(:,i) = datasample(sampMat(:,i),iters,1,'Weights',tmpWts);
    end
end
end

function [shufDistro] = bootColDs(sampMat,iters)
%Equal weight for any nonzero, non 1 value, 0 weight for 0s,1s
tmpMax = max(max(sampMat));
shufDistro = zeros(iters,size(sampMat,2));
for i = 1:size(sampMat,2)   %For number of columns
    tmpWts = sampMat(:,i) > 0;
    if sum(tmpWts) == 0   %If column is all 0s
        shufDistro(:,i) = tmpMax;
    else
    tmpWts = ones(size(sampMat,1),1).*tmpWts;   %Set wt 0 for any 1 or 0
    shufDistro(:,i) = datasample(sampMat(:,i),iters,1,'Weights',tmpWts);
    end
end
end

function [bestNLocs] = bestSeqN(nMat,pVector,nRamps)
% Looks for first index of best seq length == N
indMat = flipud(repmat(pVector,nRamps,1)');
tmpMax = max(nMat);
bestNInds = zeros(1,nRamps);
for i = 1:nRamps
    bestNInds(i) = find(nMat(:,i) == tmpMax(i),1,'last');
    if bestNInds(i) == 0
        bestNInds(i) = find(nMat(:,i) == max(nMat(:,i)),1,'first');
    end
end
bestNLocs = indMat(bestNInds);
end


