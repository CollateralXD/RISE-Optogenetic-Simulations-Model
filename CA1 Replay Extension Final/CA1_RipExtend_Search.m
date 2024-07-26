%%% Sample code for CA1 Replay Extension Model
%LKW 10/6/2021
%Relies on ripExtend_fast_V2.m and ripExtend_CohenStats.m
%Searches ripple prolongation parameter space in dimensions of ramp degree
%and input pulse duration
%Plots various metrics and statistics for evaluating network performance
%under varying waveform shape input

clearvars
close all

pStruct.rampTypeFlag        = 1;        %1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 5 = BR IMA; 7 = BR IP;
pStruct.simTypeFlag         = 2;        %1 = Linear; 2 = Linear with Adaptation
pStruct.noiseFlag           = 0;        %0 = no noise; 1 = White noise;
saveFlag                    = 0;
saveDir = 'Your\Dir\Here\';
suppGraphFlag               = 0;
disp(['Ramp Type ', num2str(pStruct.rampTypeFlag),'; Sim Type ', num2str(pStruct.simTypeFlag), '; Noise type ', num2str(pStruct.noiseFlag)]);

%Setup neuron and weight paramters
pStruct.N                   = 15;               %Nodes per region
N                           = 15;
pStruct.Ww                  = 0.0331;           %Weight strength pyr to pyr; try 0.03 for no adapt or 0.0331 for adapt
pStruct.Hh                  = 0.034;            %Weight strength IN to Pyr
pStruct.Wh                  = 0.05;             %Weight strength pyr to IN
pStruct.Wz                  = 0.02;             %Weight strength CA3 pyr to CA1 pyr
pStruct.Qz                  = 0.045;             %Weight strength CA1 IN to CA1 Pyr, 0.05 non-adapt
pStruct.Wq                  = 0.02;            %Weight strength CA3 Pyr to CA1 IN, 0.023 non-adapt
pStruct.Zq                  = 0.05;             %Weight strength CA1 Pyr to CA1 IN
pStruct.Zz                  = 0.002;            %Weight strength CA1 Pyr to CA1 Pyr
pStruct.HAuto               = 0.003;            %Weight strength IN to IN
pStruct.tha                 = 4*ones(1,N);
pStruct.thh                 = 4*ones(1,N);
pStruct.eta                 = 0.01;             %Decay constant
pStruct.actThresh           = 10*ones(1,N);     %High node firing threshold

%Setup input parameters
pStruct.T                   = 1000;             %Time steps, must be even
pStruct.cueN                = 1;
pStruct.Iexcit1             = 1;                %Cue pulse strength
pStruct.Iexcit2             = 0.1;             %Opto pulse strength
pStruct.inDur1              = 20;               %Duration for cue
pStruct.inDur2              = 100;              %Opto pulse duration
pStruct.rampLen             = 0.5*pStruct.inDur2;  %Duration ramp length e.g. 0.0 to 1.0
pStruct.onsetDelay          = 50;               %Wait time to ripple start from sim start
pStruct.stimDelay           = 130 + pStruct.inDur1 + pStruct.onsetDelay;              %Wait time to opto pulse from sim start; try 130 (+50 +20)

% Ionic Currents and related parameters 
pStruct.mu                  = 0.01;      %Ca-dependent K-current
pStruct.gm                  = 0.001;     %Gamma; voltage-dependent Ca-currents
pStruct.om                  = 0.001;     %Omega; constant for diffusion of intracellular Ca
pStruct.thc                 = 4*ones(1,N);

% Noise Related parameters
pStruct.kern                = 1;           %Kernel of random seed
pStruct.noiseAmp            = 0;         %Amplitude of Voltage noise.
pStruct.noiseMu             = 1;           %Mean of ChR2 noise distro
pStruct.noiseSigma          = 0.1;        %Variance of ChR2 noise distro

% Weight matrix
wtBias3                     = linspace(0.004,-0.00,N);
wtBias1                     = linspace(0.006,-0.006,N);
wtBias31                    = linspace(-0.01,0.01,N);
% wtBias3  = zeros(1,N);
W                           = zeros(N,N);   %Wt mat CA3 Pyr to CA3 Pyr
AH                          = zeros(N,N);   %Wt mat CA3 Pyr to CA3 IN
H                           = zeros(N,N);   %Wt mat CA3 IN to CA3 Pyr
WZ                          = zeros(N,N);   %Wt mat CA3 Pyr to CA1 Pyr
QZ                          = zeros(N,N);   %Wt mat CA1 IN to CA1 Pyr
WQ                          = zeros(N,N);   %Wt mat CA3 Pyr to CA1 IN
ZQ                          = zeros(N,N);   %Wt mat CA1 Pyr to CA1 IN
ZZ                          = ones(N,N)*pStruct.Zz;
for i = 1:N         %Build weight mats
    W(i,i)  = pStruct.Ww+wtBias3(i);   %CA3 Autorecurrency
    H(i,i)  = pStruct.Hh;   %Direct feedback IN to Pyr
    AH(i,i) = pStruct.Wh;   %Direct excitation pyr to IN
    WZ(i,i) = pStruct.Wz+wtBias1(i);   %Wt mat CA3 Pyr to CA1 Pyr
    WQ(i,i) = pStruct.Wq+wtBias31(i);  %Wt mat CA3 Pyr to CA1 IN
    QZ(i,i) = pStruct.Qz;   %Wt mat CA1 IN to CA1 Pyr
    ZQ(i,i) = pStruct.Zq;   %Wt mat CA3 Pyr to CA1 IN
    ZZ(i,i) = 0;            %CA1 Autorecurrency
    if i <= N - 1   %Forward 1
        W(i,i+1)  = (pStruct.Ww+wtBias3(i))/2;   %Pyr 2 Pyr
        WZ(i,i+1) = (pStruct.Wz+wtBias1(i))/2;
        WQ(i,i+1) = (pStruct.Wq+wtBias31(i))/2;   %Wt mat CA3 Pyr to CA1 IN
%         H(i,i+1)  = Hh/2;   %IN  2 Pyr
%         AH(i,i+1) = Wh/2;   %Pyr 2 IN
    end
    if i <= N - 2  %Forward 2
        W(i,i+2)  = (pStruct.Ww+wtBias3(i))/4;   %Pyr 2 Pyr
        WZ(i,i+2) = (pStruct.Wz+wtBias1(i))/4;
%         H(i,i+2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i+2) = Wh/4;   %Pyr 2 IN
    end
    if i > 1    %Back 1
        WQ(i,i-1) = (pStruct.Wq+wtBias31(i))/2;   %Wt mat CA3 Pyr to CA1 IN
%         W(i,i-1)  = (Ww+wtBias(i))/2;   %Pyr 2 Pyr
%         H(i,i-1)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-1) = Wh/2;   %Pyr 2 IN
    end
    if i > 2    %Back 2   
%         W(i,i-2)  = (Ww+wtBias(i))/4;   %Pyr 2 Pyr
%         H(i,i-2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-2) = Wh/4;   %Pyr 2 IN
    end
end

pStruct.W = W; pStruct.H = H; pStruct.AH = AH;  pStruct.ZZ = ZZ; %Add wtmats to pStruct
pStruct.WZ = WZ; pStruct.QZ = QZ; pStruct.WQ = WQ; pStruct.ZQ = ZQ;

%Set up search variables and parameters
nPs                     = 51;
pLim                    = 250;      %250 for inDur2; 500 for stimDelay; 1 for noiseSigma
pFloor                  = 0;
nRamps                  = 21;
pVect                   = linspace(pFloor,pLim,nPs);    %Vector of parameter values
if mod(pStruct.rampTypeFlag,2) == 1; rampLim = 1; else; rampLim = 0.5; end %Set ramp limit to 1 or 0.5 depending on flag
rampPercs = linspace(0,rampLim,nRamps);

outputs.ds          = zeros(nRamps,nPs);
outputs.cs          = cell(nRamps,nPs);
outputs.means       = zeros(nRamps,nPs);
outputs.SDs         = zeros(nRamps,nPs);
outputs.Ns          = zeros(nRamps,nPs);

%% Gut check above parameters
[actCell,hactCell,inCell] = ripEx_CA1(pStruct);    %Using ramp defined above
ca3P = actCell{1}; ca1PRmp = actCell{2}; ca1PCtl = actCell{3};
ca3I = hactCell{1}; ca1IRmp = hactCell{2}; ca1ICtl = hactCell{3};
ca3A = inCell{1}; ca1ARmp = inCell{2}; ca1ACtl = inCell{3};
dGut = ripEx_CA1_cohen(actCell,hactCell,pStruct);

% Calculate Peaks
locsR = []; locsC = [];
actThresh = 10;
for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(ca1PRmp(i,:));
    if ~isempty(locsTmp) && pksTmp(1) > actThresh
        indtmp = find(ca1PRmp(i,:)>actThresh,1);
        locsR = [locsR indtmp]; 
    end
    [pksTmp,locsTmp] = findpeaks(ca1PCtl(i,:));
    if ~isempty(locsTmp) && pksTmp(1) > actThresh
        indtmp = find(ca1PCtl(i,:)>actThresh,1);
        locsC = [locsC indtmp]; 
    end
end

gutRmpPlot = plotCA1RipEx(ca3P,ca1PRmp,ca1ARmp,locsR,pStruct);
gutCtlPlot = plotCA1RipEx(ca3P,ca1PCtl,ca1ACtl,locsC,pStruct);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.5, 0.5, 0.38]);
legend('CA3 Pyr','CA1 Ctl Pyr')

%% Run Block
for i = 1:nRamps
    pStruct.tmpRamp = rampPercs(i);         %Define ramp percentage
    for j = 1:nPs
        pStruct.inDur2      = pVect(j); 	  %2nd search param: inDur
%         pStruct.noiseSigma  = pVect(j);     %2nd search param: noiseSigma of ChR2 expression
%         pStruct.noiseAmp    = pVect(j);     %2nd search param: amplitude of voltage noise
%         pStruct.stimDelay   = pVect(j) + pStruct.onsetDelay + pStruct.inDur1; %2nd search param: stimDelay
        pStruct.rampLen     = round(pStruct.tmpRamp*pStruct.inDur2);   %Define ramp length
        %Run function
        [actTmp,hactTmp,inTmp] = ripEx_CA1(pStruct);     %Run simulation
        %Calculate stats
        [outputs.ds(i,j),~,outputs.dINs(i,j),outputs.Ns(i,j)] = ripEx_CA1_cohen(actTmp,hactTmp,pStruct);
    end
end

%% Clean data
outputs.dsNanCor = abs(outputs.ds);
if pStruct.noiseFlag == 1   %Use 99th percentile instead of data max
    linearDs = abs(reshape(outputs.ds,[1 numel(outputs.ds)]));
    tmpMax = prctile(linearDs,99)
else                        %Use maximum, non Inf simulation value
    tmpMax = max(max(abs(~isinf(outputs.dsNanCor).*outputs.dsNanCor)));
end
outputs.dsNanCor(isnan(outputs.dsNanCor)) = tmpMax;
outputs.dsNanCor(isinf(outputs.dsNanCor)) = tmpMax; 
outputs.dINs = abs(outputs.dINs); tmpINMax = max(max(outputs.dINs)); 
outputs.dINs(isnan(outputs.dINs)) = tmpINMax;

% Other metrics
minDLocs = minTimes(outputs.dsNanCor,pVect,nRamps); %Calculate duration of minimum effect size
minDs = min(flipud(outputs.dsNanCor'));    %Calculate minimum effect size

% cbMax = tmpMax;
% cbMax = max(max(outputs.dsNanCor));
cbMax = 2.5;
%% Plotting of Heatmaps

%Plot map of d-scores
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaa = gca;
imagesc(flipud(outputs.dsNanCor'),[0,cbMax]);
aaa.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaa.XTick = 1:round(nRamps-1)/5:nRamps;
yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
colormap hot; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Opto vs Control");
% else; title("Double Ramp Opto vs Control"); end
% hcb.Label.String = "Cohen's d";
set(aaa,'FontSize',24,'fontname','times')

%Plot map of sequence lengths
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaz = gca;
imagesc(flipud(abs(outputs.Ns')),[0,N]);
aaz.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaz.XTick = 1:round(nRamps-1)/5:nRamps;
yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
axis square; hcb = colorbar;
% hcb.Label.String = "Sequence Length";
set(aaz,'FontSize',24,'fontname','times')

% % Plot map of IN d-scores
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aax = gca;
% imagesc(flipud(outputs.dINs'),[0,max(max(outputs.dINs))]);
% aax.YTick = 1:round((nPs-1)/5):nPs; 
% aax.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% % xlabel('Ramp Percentage'); 
% % ylabel('Learn Overlap (ms)'); 
% colormap cool; axis square; hcb = colorbar;
% set(aax,'FontSize',24,'fontname','times')

%% Shuffles
permN = 1000;
z = 1.96;   %confidence level 95% = 1.96, 99% = 2.58
durShufDistro = datasample(flipud(outputs.dsNanCor'),permN,1);
INShufDistro  = datasample(flipud(outputs.dINs'),permN,1);
seqShufDistro = datasample(flipud(outputs.Ns'),permN,1);
rmpShufDistro = datasample(flipud(outputs.dsNanCor),permN,1);
rmpSeqShufDistro = datasample(flipud(outputs.Ns),permN,1);
shuf.durShuf = mean(durShufDistro);
shuf.INShuf = mean(INShufDistro);
shuf.seqShuf = mean(seqShufDistro);
shuf.rmpShuf = mean(rmpShufDistro);
shuf.rmpSeqShuf = mean(rmpSeqShufDistro);
shuf.sdDurShuffle = std(durShufDistro); 
shuf.sdINShuffle  = std(INShufDistro);
shuf.sdSeqShuffle = std(seqShufDistro);
shuf.sdRmpShuffle = std(rmpShufDistro);
shuf.sdRmpSeqShuffle = std(rmpSeqShufDistro);

% 95% CIs
shuf.upCIdur = shuf.durShuf + z*shuf.sdDurShuffle/sqrt(permN);
shuf.dnCIdur = shuf.durShuf - z*shuf.sdDurShuffle/sqrt(permN);
shuf.upCIseq = shuf.seqShuf + z*shuf.sdSeqShuffle/sqrt(permN);
shuf.dnCIseq = shuf.seqShuf - z*shuf.sdSeqShuffle/sqrt(permN);
shuf.upCIIN  = shuf.INShuf  + z*shuf.sdINShuffle/sqrt(permN);
shuf.dnCIIN  = shuf.INShuf  - z*shuf.sdINShuffle/sqrt(permN);
shuf.upCIrmp = shuf.rmpShuf + z*shuf.sdRmpShuffle/sqrt(permN);
shuf.dnCIrmp = shuf.rmpShuf - z*shuf.sdRmpShuffle/sqrt(permN);
shuf.upCIrmpSeq = shuf.rmpSeqShuf + z*shuf.sdRmpSeqShuffle/sqrt(permN);
shuf.dnCIrmpSeq = shuf.rmpSeqShuf - z*shuf.sdRmpSeqShuffle/sqrt(permN);

% Plot Basic Shuffle Figures
rVect = rampPercs*100;
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); ccb = gca;
hold on;
plot(rVect,shuf.durShuf,'k',rVect,shuf.upCIdur,'k--',rVect,shuf.dnCIdur,'k--')
% plot(rVect,shuf.INShuf,'r',rVect,shuf.upCIIN,'r--',rVect,shuf.dnCIIN,'r--')
% plot(pVect,shuf.rmpShuf,'b',pVect,shuf.upCIrmp,'b--',pVect,shuf.dnCIrmp,'b--')
% plot(pVect,shuf.rmpSeqShuf,'c',pVect,shuf.upCIrmpSeq,'c--',pVect,shuf.dnCIrmpSeq,'c--')
% plot(rVect,shuf.seqShuf,'b',rVect,shuf.upCIseq,'b--',rVect,shuf.dnCIseq,'b--')
if pStruct.rampTypeFlag == 1; legend('Pyr Shuffle','upper CI','lower CI','IN Shuffle'); end
ylim([0 cbMax]);
set(gca,'FontSize',24,'fontname','times')
% ylabel('Mean Effect Size'); xlabel('Ramp Percentage')
ylabel('Mean Effect Size'); xlabel('2nd search parameter')

% % Minimum effect sizes
% figure;
% plot(rampPercs*100,minDs)

%% Supplementary Graphs
set(0,'DefaultLineLineWidth',2)

if suppGraphFlag == 1
    
% % Afferent Input 3D plot code
% sqGhost = zeros(1,pStruct.T);
% sqGhost(pStruct.stimDelay+1:pStruct.stimDelay+pStruct.inDur2) = pStruct.Iexcit2;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.35, 0.45]); gaa = gca;
% w3 = waterfall(1:pStruct.T,1:N,ARamp); hold on;
% w3.LineWidth = 3;
% % w3Sq = waterfall(1:pStruct.T,1,sqGhost);
% % w3Sq.LineStyle = '--'; w3Sq.EdgeColor = 'r'; w3Sq.LineWidth = 2;
% % xlabel('Time')
% gaa.YTick = 1:7:15;
% % ylabel('Node'); zlabel('Amplitude')
% set(gca,'FontSize',24,'fontname','times')

% Combined pulse duration Shuffles - must load in previous data
patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); saa = gca;
axis square; hold on;
plot(FR_IMA_shufs.durShuf,'r-');
plot(DR_IMA_shufs.durShuf,'b-'); 
plot(BR_IMA_shufs.durShuf,'c-');
plot(FR_IP_shufs.durShuf,'r--');
plot(DR_IP_shufs.durShuf,'b--');
plot(BR_IP_shufs.durShuf,'c--');
patch(patchXs,[FR_IMA_shufs.dnCIdur,fliplr(FR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[DR_IMA_shufs.dnCIdur,fliplr(DR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IMA_shufs.dnCIdur,fliplr(BR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[FR_IP_shufs.dnCIdur,fliplr(FR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[DR_IP_shufs.dnCIdur,fliplr(DR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IP_shufs.dnCIdur,fliplr(BR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
ylim([0 3]); xlim([1 nRamps]); saa.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
set(gca,'FontSize',24,'fontname','times')
% 
% % Combined IN vs Effect Size Shuffles - must load in previous data
% patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sab = gca;
% axis square; hold on;
% plot(FR_IMA_shufs.INShuf,'r-');
% plot(DR_IMA_shufs.INShuf,'b-'); 
% plot(BR_IMA_shufs.INShuf,'c-');
% plot(FR_IP_shufs.INShuf,'r--');
% plot(DR_IP_shufs.INShuf,'b--');
% plot(BR_IP_shufs.INShuf,'c--');
% patch(patchXs,[FR_IMA_shufs.dnCIIN,fliplr(FR_IMA_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[DR_IMA_shufs.dnCIIN,fliplr(DR_IMA_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IMA_shufs.dnCIIN,fliplr(BR_IMA_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[FR_IP_shufs.dnCIIN,fliplr(FR_IP_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[DR_IP_shufs.dnCIIN,fliplr(DR_IP_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IP_shufs.dnCIIN,fliplr(BR_IP_shufs.upCIIN)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','southwest');
% ylim([0 3]); xlim([1 nRamps]); sab.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')

% Combined SeqLen Shuffles - must load in previous data
patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sac = gca;
axis square; hold on;
plot(FR_IMA_shufs.seqShuf,'r-');
plot(DR_IMA_shufs.seqShuf,'b-'); 
plot(BR_IMA_shufs.seqShuf,'c-');
plot(FR_IP_shufs.seqShuf,'r--');
plot(DR_IP_shufs.seqShuf,'b--');
plot(BR_IP_shufs.seqShuf,'c--');
patch(patchXs,[FR_IMA_shufs.dnCIseq,fliplr(FR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[DR_IMA_shufs.dnCIseq,fliplr(DR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IMA_shufs.dnCIseq,fliplr(BR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[FR_IP_shufs.dnCIseq,fliplr(FR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[DR_IP_shufs.dnCIseq,fliplr(DR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IP_shufs.dnCIseq,fliplr(BR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
xlabel('Ramp Percentage'); ylabel('Mean Sequence Length')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','southwest');
ylim([5 15]); xlim([1 nRamps]); sac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
set(gca,'FontSize',24,'fontname','times')

% Combined Ramp Percentage Shuffles: Perturbation vs Delay Length - must load in previous data
patchXs = [1:nPs,linspace(nPs,1,nPs)]; %Vector of x coords;
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sad = gca;
axis square; hold on;
plot(FR_IMA_shufs.rmpShuf,'r-');
plot(DR_IMA_shufs.rmpShuf,'b-'); 
plot(BR_IMA_shufs.rmpShuf,'c-');
plot(FR_IP_shufs.rmpShuf,'r--');
plot(DR_IP_shufs.rmpShuf,'b--');
plot(BR_IP_shufs.rmpShuf,'c--');
patch(patchXs,[FR_IMA_shufs.dnCIrmp,fliplr(FR_IMA_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[DR_IMA_shufs.dnCIrmp,fliplr(DR_IMA_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IMA_shufs.dnCIrmp,fliplr(BR_IMA_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[FR_IP_shufs.dnCIrmp,fliplr(FR_IP_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[DR_IP_shufs.dnCIrmp,fliplr(DR_IP_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IP_shufs.dnCIrmp,fliplr(BR_IP_shufs.upCIrmp)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
xlabel('Input Duration'); ylabel('Mean Effect Size')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','southeast');
ylim([0 inf]); xlim([1 nPs]); sad.XTick = linspace(0,nPs,6); xticklabels(linspace(0,pLim,6))
set(gca,'FontSize',24,'fontname','times')

% Combined Ramp Percentage Shuffles: Sequence Length vs Delay Length - must load in previous data
patchXs = [1:nPs,linspace(nPs,1,nPs)]; %Vector of x coords;
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); sae = gca;
axis square; hold on;
plot(FR_IMA_shufs.rmpSeqShuf,'r-');
plot(DR_IMA_shufs.rmpSeqShuf,'b-'); 
plot(BR_IMA_shufs.rmpSeqShuf,'c-');
plot(FR_IP_shufs.rmpSeqShuf,'r--');
plot(DR_IP_shufs.rmpSeqShuf,'b--');
plot(BR_IP_shufs.rmpSeqShuf,'c--');
patch(patchXs,[FR_IMA_shufs.dnCIrmpSeq,fliplr(FR_IMA_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[DR_IMA_shufs.dnCIrmpSeq,fliplr(DR_IMA_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IMA_shufs.dnCIrmpSeq,fliplr(BR_IMA_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[FR_IP_shufs.dnCIrmpSeq,fliplr(FR_IP_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[DR_IP_shufs.dnCIrmpSeq,fliplr(DR_IP_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
patch(patchXs,[BR_IP_shufs.dnCIrmpSeq,fliplr(BR_IP_shufs.upCIrmpSeq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
xlabel('Input Duration'); ylabel('Mean Sequence Length')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northeast');
ylim([0 15]); xlim([1 nPs]); sae.XTick = linspace(0,nPs,6); xticklabels(linspace(0,pLim,6))
set(gca,'FontSize',24,'fontname','times')

%Must Load in previously saved minDLoc data
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); saf = gca;
plot(FR_IMA_minLocs,'r-'); hold on;
plot(DR_IMA_minLocs,'b-');
plot(BR_IMA_minLocs,'c-');
plot(FR_IP_minLocs,'r--');
plot(DR_IP_minLocs-1,'b--');
plot(BR_IP_minLocs+1,'c--');
saf.XTick = 1:round(nRamps-1)/5:nRamps;
xticklabels(0:20:100);
xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
axis square; 
% title('Input Duration of Minimum Effect Size')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
ylim([0 101]); xlim([1 nRamps]);
set(saf,'FontSize',24,'fontname','times')

%Must Load in previously saved minD data
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); sag = gca;
plot(FR_IMA_minDs,'r-'); hold on;
plot(DR_IMA_minDs,'b-');
plot(BR_IMA_minDs,'c-');
plot(FR_IP_minDs,'r--');
plot(DR_IP_minDs,'b--');
plot(BR_IP_minDs,'c--');
sag.XTick = 1:round(nRamps-1)/5:nRamps;
xticklabels(0:20:100);
xlabel('Ramp Percentage'); ylabel('Least Temporal Disruption'); 
axis square; 
% title('Minimum Effect Size')
legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
ylim([0 0.5]); xlim([1 nRamps]);
set(sag,'FontSize',24,'fontname','times')

    if saveFlag == 1
        saveDir = 'F:\Research\Code\CA1 Region Code\CA1_RipExtend\CA1_RipExtend_Figs\';
        saveas(saa,[saveDir,'allShufs_rampsXmeanDs',figtag],'fig')
        saveas(saa,[saveDir,'allShufs_rampsXmeanDs',figtag],'png')
        saveas(saa,[saveDir,'allShufs_rampsXmeanDs',figtag],'svg')
%         saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'fig')
%         saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'png')
%         saveas(sab,[saveDir,'allShufs_rampsXminDs',figtag],'svg')
        saveas(sac,[saveDir,'allShufs_rampsXmeanNs',figtag],'fig')
        saveas(sac,[saveDir,'allShufs_rampsXmeanNs',figtag],'png')
        saveas(sac,[saveDir,'allShufs_rampsXmeanNs',figtag],'svg')
        saveas(sad,[saveDir,'allShufs_inDursXmeanDs',figtag],'fig')
        saveas(sad,[saveDir,'allShufs_inDursXmeanDs',figtag],'png')
        saveas(sad,[saveDir,'allShufs_inDursXmeanDs',figtag],'svg')
        saveas(sae,[saveDir,'allShufs_inDursXmeanNs',figtag],'fig')
        saveas(sae,[saveDir,'allShufs_inDursXmeanNs',figtag],'png')
        saveas(sae,[saveDir,'allShufs_inDursXmeanNs',figtag],'svg')
        saveas(saf,[saveDir,'allShufs_rampsXminLocs',figtag],'fig')
        saveas(saf,[saveDir,'allShufs_rampsXminLocs',figtag],'png')
        saveas(saf,[saveDir,'allShufs_rampsXminLocs',figtag],'svg')        
        saveas(sag,[saveDir,'allShufs_rampsXminDs',figtag],'fig')
        saveas(sag,[saveDir,'allShufs_rampsXminDs',figtag],'png')
        saveas(sag,[saveDir,'allShufs_rampsXminDs',figtag],'svg')        
    end
    
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

end
%% Saves
if saveFlag == 1
disp('saving vars and figs')

if z == 1.96
    ciStr = ['_perms',num2str(permN),'_CI95'];
elseif z == 2.58
    ciStr = ['_perms',num2str(permN),'_CI99'];
end

prmtag = '_Iexcit2_0_';            %name of fixed parameter
prmtmp = 100*pStruct.Iexcit2;      %value of fixed parameter
figtag = [prmtag,num2str(prmtmp)]; %name-value string for figure saves
% figtag = ['_adapt_Ww_',1000*num2str(pStruct.Ww)];

if pStruct.rampTypeFlag == 1
    sBase = ['FR_IMA',figtag,ciStr];
    FR_IMA_minLocs  = minDLocs;
    FR_IMA_minDs    = minDs;
    FR_IMA_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName,'FR_IMA_minLocs','FR_IMA_minDs','FR_IMA_shufs')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'FR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 2
    sBase = ['DR_IMA',figtag,ciStr];
    DR_IMA_minLocs  = minDLocs;
    DR_IMA_minDs    = minDs;
    DR_IMA_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName, 'DR_IMA_minLocs','DR_IMA_minDs','DR_IMA_shufs')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'DR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 5
    sBase = ['BR_IMA',figtag,ciStr];
    BR_IMA_minLocs  = minDLocs;
    BR_IMA_minDs    = minDs;
    BR_IMA_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName, 'BR_IMA_minLocs','BR_IMA_minDs','BR_IMA_shufs')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'BR_IMA_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 3
    sBase = ['FR_IP',figtag,ciStr];
    FR_IP_minLocs  = minDLocs;
    FR_IP_minDs    = minDs;
    FR_IP_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName, 'FR_IP_minLocs','FR_IP_minDs','FR_IP_shufs')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'FR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'FR_IP_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 4
    sBase = ['DR_IP',figtag,ciStr];
    DR_IP_minLocs  = minDLocs;
    DR_IP_minDs    = minDs;
    DR_IP_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName, 'DR_IP_minLocs','DR_IP_minDs','DR_IP_shufs')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'DR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'DR_IP_heatmap_seqs',figtag],'svg')
elseif pStruct.rampTypeFlag == 7
    sBase = ['BR_IP',figtag,ciStr];
    BR_IP_minLocs  = minDLocs;
    BR_IP_minDs    = minDs;
    BR_IP_shufs    = shuf;
    fName = [saveDir, sBase];
    save(fName, 'BR_IP_minLocs','BR_IP_minDs','BR_IP_shufs')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'fig')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'png')
    saveas(aaa,[saveDir,'BR_IP_heatmap_ds',figtag],'svg')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'fig')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'png')
    saveas(aaz,[saveDir,'BR_IP_heatmap_seqs',figtag],'svg')
end
disp('Saves done.')
end
%% Functions

function [minLocs] = minTimes(dsMat,pVector,nRamps)
dsMat = flipud(abs(dsMat'));
[~,dInds] = min(dsMat,[],1);
%     minBin = dsMat == dMins;
indMat = flipud(repmat(pVector,nRamps,1)');
%     minCol = minBin.*indMat;     %Removes all but earliest minimum
minLocs = indMat(dInds);
end



