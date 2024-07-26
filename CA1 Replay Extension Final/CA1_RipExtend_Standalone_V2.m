%%% Sample code for Replay Extension Model
%LKW 10/1/2021
%Relies on plotRipExtendSingle.m for plots
%Calculates activity of CA3 pyr circuit undergoing a cued ripple/replay
%event with or without a secondary optogenetic pulse applied mid-replay. 
%Incorporates simple adaptation current using simulated intracellular Calcium dynamics
%Plots rate-based circuit activity. 
%%%

clearvars
close all
% rng(3)                          %For reproducibility of stochastic simulations

rampTypeFlag        = 1;        %'ctl' = no 2nd pulse; 1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 5 = BR IMA; 7 = BR = IP
simTypeFlag         = 2;        %1 = Linear; 2 = Linear with Adaptation
noiseFlag           = 0;        %0 = no noise; 1 = White noise
disp(['Ramp Type ', num2str(rampTypeFlag),'; Sim Type ', num2str(simTypeFlag), '; Noise type ', num2str(noiseFlag)]);
saveFlag            = 0;
saveDir = 'Your\Dir\Here\';

cueN                = 1;            %Cue node
N                   = 15;           %Nodes per region
Ww                  = 0.0331;       %Weight strength CA3 pyr to CA3 pyr, try 0.03 for non-adapt, 0.0331 for adapt
Hh                  = 0.034;        %Weight strength CA3 IN to CA3 Pyr, try 0.034 non-adapt
Wh                  = 0.05;         %Weight strength CA3 pyr to CA3 IN, try 0.05 non-adapt
Wz                  = 0.02;         %Weight strength CA3 pyr to CA1 pyr, try 0.02 non-adapt
Qz                  = 0.045;        %Weight strength CA1 IN to CA1 Pyr, try 0.05 non-adapt, 0.045 for adapt
Wq                  = 0.02;         %Weight strength CA3 Pyr to CA1 IN, try 0.023 non-adapt, 0.02 for adapt
Zq                  = 0.05;         %Weight strength CA1 Pyr to CA1 IN, try 0.05 non-adapt
Zz                  = 0.002;        %Weight strength CA1 Pyr to CA1 Pyr, try 0.002 non-adapt
HAuto               = 0.003;        %Inhibition self feedback, try 0.003
tha                 = 4*ones(1,N);
thh                 = 4*ones(1,N);
thc                 = 4*ones(1,N);
eta                 = 0.01;         %Decay constant
T                   = 1500;         %Time steps, must be even
Iexcit1             = 1;            %Try 0.5 with Iexcit2 at 0.05 for ramp < square ripple
Iexcit2             = 0.1;         %Try as low as 0.09 or 0.05 opto pulse for good ripple spacing
inDur1              = 20;           %Duration for kicking off a ripple
inDur2              = 100;          %Stim duration
rampPerc            = 0.25;         %Percentage of ramp
rampLen             = round(inDur2*rampPerc);  %Duration ramp length e.g. 1/5 or 1/2
onsetDelay          = 50;           %Wait time to ripple start from sim start
stimDelay           = 130 + onsetDelay + inDur1;          %Wait time to opto pulse from simulation start
noiseAmp            = 0;         %Amplitude of Voltage noise.
noiseMu             = 1;           %Mean of ChR2 noise distro
noiseSigma          = 0.05;         %Variance of ChR2 noise distro

% Ionic Currents and related parameters 
Ek  = -10;              %Reversal Potential of Potassium
Ena = 70;               %Sodium
Ecl = 0;                %Chlorine
mu  = 0.01;             %Ca-dependent K-current
gm  = 0.001;            %Gamma; voltage-dependent Ca-currents
om  = 0.001;            %Omega; constant for diffusion of intracellular Ca

% Weight matrix
wtBias3  = linspace(0.004,-0.00,N);
wtBias1  = linspace(0.006,-0.006,N);
wtBias31 = linspace(-0.01,0.01,N);
% wtBias3  = zeros(1,N);
W       = zeros(N,N);   %Wt mat CA3 Pyr to CA3 Pyr
AH      = zeros(N,N);   %Wt mat CA3 Pyr to CA3 IN
H       = zeros(N,N);   %Wt mat CA3 IN to CA3 Pyr
WZ      = zeros(N,N);   %Wt mat CA3 Pyr to CA1 Pyr
QZ      = zeros(N,N);   %Wt mat CA1 IN to CA1 Pyr
WQ      = zeros(N,N);   %Wt mat CA3 Pyr to CA1 IN
ZQ      = zeros(N,N);   %Wt mat CA1 Pyr to CA1 IN
ZZ      = ones(N,N)*Zz;   %Wt mat CA1 Pyr to CA1 Pyr

for i = 1:N         %Build weight mats
    W(i,i)  = Ww+wtBias3(i);   %CA3 Autorecurrency
    H(i,i)  = Hh;   %Direct feedback IN to Pyr
    AH(i,i) = Wh;   %Direct excitation pyr to IN
    WZ(i,i) = Wz+wtBias1(i);   %Wt mat CA3 Pyr to CA1 Pyr
    WQ(i,i) = Wq+wtBias31(i);   %Wt mat CA3 Pyr to CA1 IN
    QZ(i,i) = Qz; %Wt mat CA1 IN to CA1 Pyr
    ZQ(i,i) = Zq;   %Wt mat CA1 Pyr to CA1 IN
    ZZ(i,i) = 0;    %CA1 Autorecurrency
    if i <= N - 1   %Forward 1
        W(i,i+1)  = (Ww+wtBias3(i))/2;   %Pyr 2 Pyr
        WZ(i,i+1) = (Wz+wtBias1(i))/2;
        WQ(i,i+1) = (Wq+wtBias31(i))/2;   %Wt mat CA3 Pyr to CA1 IN
%         ZZ(i,i+1) = Zz;
%         H(i,i+1)  = Hh/2;   %IN  2 Pyr
%         AH(i,i+1) = Wh/2;   %Pyr 2 IN
    end
    if i <= N - 2  %Forward 2
        W(i,i+2)  = (Ww+wtBias3(i))/4;   %Pyr 2 Pyr
        WZ(i,i+2) = (Wz+wtBias1(i))/4;
%         ZZ(i,i+2) = Zz;
%         H(i,i+2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i+2) = Wh/4;   %Pyr 2 IN
    end
    if i > 1    %Back 1
        WQ(i,i-1) = (Wq+wtBias31(i))/2;   %Wt mat CA3 Pyr to CA1 IN
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

%CA3 Region
ca3A    = zeros(N,T);
ca3P    = zeros(N,T);
ca3I    = zeros(N,T);
ca3C    = zeros(N,T);
%CA1 Ramp vectors
ca1ARmp = zeros(N,T);
ca1PRmp = zeros(N,T);
ca1IRmp = zeros(N,T);
ca1CRmp = zeros(N,T);
%CA1 Square vectors
ca1ASqr = zeros(N,T);
ca1PSqr = zeros(N,T);
ca1ISqr = zeros(N,T);
ca1CSqr = zeros(N,T);
%CA1 Control vectors
ca1ACtl = zeros(N,T);
ca1PCtl = zeros(N,T);
ca1ICtl = zeros(N,T);
ca1CCtl = zeros(N,T);

%Set up waveform Afferent Inputs
ca3A(cueN,onsetDelay+1:onsetDelay+inDur1)       = Iexcit1; %Start CA3 ripple

ca1ASqr(:,stimDelay+1:stimDelay+inDur2) = Iexcit2; % rand(N,inDur);  %Square pulse / random

if rampTypeFlag == 1        %FR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen)= repmat(linspace(0,Iexcit2,rampLen),[N,1]);   %Front Ramp
elseif rampTypeFlag == 2    %DR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2; 
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,Iexcit2,rampLen),[N,1]);  %Front Ramp
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
elseif rampTypeFlag == 5    %BR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
end
squarea = Iexcit2*inDur2;     %Area under the square pulse
if rampTypeFlag == 3        %FR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
elseif rampTypeFlag == 4    %DR IP
    IRamp = squarea/(inDur2 - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
elseif rampTypeFlag == 7    %BR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
end

%% Build noise
if noiseFlag == 0
    ANoise = zeros(N,T);
    tissNoise = ones(N,1);
elseif noiseFlag == 1
%     %Voltage Noise from LFP fluctuations
    ANoise = rand(N,T)*noiseAmp - 0.5*noiseAmp;
    
    %Light Scattering Noise from distance in 3D space
%     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
%     load ScatteringFitIrr8mW.mat    %fspline calculated from 8mW irradiance (mW/mm^2)
    load ScatteringFitIrr10mW.mat    %fspline calculated from 10mW irradiance (mW/mm^2)
    sTip = [0 0 0]; %set laser source at origin
    unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
    unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
    distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
    distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
    distNoise = distNoise/5;        %Get normalized % activation with distance dropoff
    ChR2Noise = ones(N,1);
    
%     %Protein Expression Noise
%     ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression
%     distNoise = ones(N,1);
    
    tissNoise = ChR2Noise .* distNoise;   %Multiply ChR2 Noise with distance Noise
end

% %Plot histogram of unit distance from source in 3D space
% figure; 
% h11 = histogram(unitDist);
% h11.Normalization = 'probability';
% xlabel('Distance (mm) from Source'); ylabel('Probability')
% set(gca,'FontSize',20,'fontname','times')

% %Plot histogram of power dropoff from source in 3D space
% figure; 
% hold on; 
% h2 = histogram(distNoise_Irr10);
% h2.BinWidth = 0.01;
% h2.Normalization = 'probability';
% h1 = histogram(distNoise_Irr8);
% h1.BinWidth = 0.01;
% h1.Normalization = 'probability';
% xlabel('% Irradiance (mW/mm^2) of 5mW'); ylabel('Probability')
% set(gca,'FontSize',20,'fontname','times')
% 
% %Plot histogram of ChR2 expression variability
% figure; hold on;
% h22 = histogram(ChR2Noise_sigma0_05);
% h22.Normalization = 'probability';
% h22.BinWidth = 0.01;
% h23 = histogram(ChR2Noise_sigma0_1);
% h23.Normalization = 'probability';
% h23.BinWidth = 0.01;
% xlabel('% Opsin Expression Efficiency'); ylabel('Probability')
% set(gca,'FontSize',20,'fontname','times')

%% Run Block
for t=1:T-1
    for j = 1:N
        %CA3 Region
        if (ca3P(j,t)    - tha(j)) > 0; aa3(j)    = ca3P(j,t)    - tha(j); else aa3(j)    = 0; end %Threshold activity over 'tha' to 0
        if (ca3I(j,t)    - thh(j)) > 0; hh3(j)    = ca3I(j,t)    - thh(j); else hh3(j)    = 0; end %Threshold activity over 'thh' to 0
        if (ca3P(j,t)    - thc(j)) > 0; cc3(j)    = ca3P(j,t)    - thc(j); else cc3(j)    = 0; end %Threshold activity under 'thc' to 0
        %CA1 Region
        if (ca1PRmp(j,t) - tha(j)) > 0; aa1Rmp(j) = ca1PRmp(j,t) - tha(j); else aa1Rmp(j) = 0; end %Threshold activity over 'tha' to 0
        if (ca1PSqr(j,t) - tha(j)) > 0; aa1Sqr(j) = ca1PSqr(j,t) - tha(j); else aa1Sqr(j) = 0; end
        if (ca1PCtl(j,t) - tha(j)) > 0; aa1Ctl(j) = ca1PCtl(j,t) - tha(j); else aa1Ctl(j) = 0; end
        if (ca1IRmp(j,t) - thh(j)) > 0; hh1Rmp(j) = ca1IRmp(j,t) - thh(j); else hh1Rmp(j) = 0; end %Threshold activity over 'thh' to 0
        if (ca1ISqr(j,t) - thh(j)) > 0; hh1Sqr(j) = ca1ISqr(j,t) - thh(j); else hh1Sqr(j) = 0; end
        if (ca1ICtl(j,t) - thh(j)) > 0; hh1Ctl(j) = ca1ICtl(j,t) - thh(j); else hh1Ctl(j) = 0; end
        if (ca1PRmp(j,t) - thc(j)) > 0; cc1Rmp(j) = ca1PRmp(j,t) - thc(j); else cc1Rmp(j) = 0; end %Threshold activity under 'thc' to 0
        if (ca1PSqr(j,t) - thc(j)) > 0; cc1Sqr(j) = ca1PSqr(j,t) - thc(j); else cc1Sqr(j) = 0; end
        if (ca1PCtl(j,t) - thc(j)) > 0; cc1Ctl(j) = ca1PCtl(j,t) - thc(j); else cc1Ctl(j) = 0; end
    end
    
    if simTypeFlag == 1
        %Standard Linear without Adapatation
        %CA3
        da3    = tissNoise.*ca3A(:,t) + ANoise(:,t) + (aa3*W)' - (hh3*H)' - eta.*ca3P(:,t);
        dh3    = (aa3*AH)' - HAuto.*hh3' - eta.*ca3I(:,t);
        %CA1
        da1Rmp = tissNoise.*ca1ARmp(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Rmp*ZZ)' - (hh1Rmp*QZ)' - eta.*ca1PRmp(:,t);
        dh1Rmp = (aa3*WQ)' + (aa1Rmp*ZQ)' - HAuto.*hh1Rmp' - eta.*ca1IRmp(:,t);
        da1Sqr = tissNoise.*ca1ASqr(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Sqr*ZZ)' - (hh1Sqr*QZ)' - eta.*ca1PSqr(:,t);
        dh1Sqr = (aa3*WQ)' + (aa1Sqr*ZQ)' - HAuto.*hh1Sqr' - eta.*ca1ISqr(:,t);
        da1Ctl = tissNoise.*ca1ACtl(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Ctl*ZZ)' - (hh1Ctl*QZ)' - eta.*ca1PCtl(:,t);
        dh1Ctl = (aa3*WQ)' + (aa1Ctl*ZQ)' - HAuto.*hh1Ctl' - eta.*ca1ICtl(:,t);
    elseif simTypeFlag == 2
        da3    = tissNoise.*ca3A(:,t) + ANoise(:,t) + (aa3*W)' - (hh3*H)' - eta.*ca3P(:,t) + mu.*ca3C(:,t).*(Ek - ca3P(:,t));
        dh3    = (aa3*AH)' - HAuto.*hh3' - eta.*ca3I(:,t);
        dc3    = gm.*cc3' - om.*ca3C(:,t);
        %Linear variant with Adaptation
        da1Rmp = tissNoise.*ca1ARmp(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Rmp*ZZ)' - (hh1Rmp*QZ)' - eta.*ca1PRmp(:,t) + mu.*ca1CRmp(:,t).*(Ek - ca1PRmp(:,t));
        dh1Rmp = (aa3*WQ)' + (aa1Rmp*ZQ)' - HAuto.*hh1Rmp' - eta.*ca1IRmp(:,t);
        dc1Rmp = gm.*cc1Rmp' - om.*ca1CRmp(:,t);
        da1Sqr = tissNoise.*ca1ASqr(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Sqr*ZZ)' - (hh1Sqr*QZ)' - eta.*ca1PSqr(:,t) + mu.*ca1CSqr(:,t).*(Ek - ca1PSqr(:,t));
        dh1Sqr = (aa3*WQ)' + (aa1Sqr*ZQ)' - HAuto.*hh1Sqr' - eta.*ca1ISqr(:,t);
        dc1Sqr = gm.*cc1Sqr' - om.*ca1CSqr(:,t);
        da1Ctl = tissNoise.*ca1ACtl(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Ctl*ZZ)' - (hh1Ctl*QZ)' - eta.*ca1PCtl(:,t) + mu.*ca1CCtl(:,t).*(Ek - ca1PCtl(:,t));
        dh1Ctl = (aa3*WQ)' + (aa1Ctl*ZQ)' - HAuto.*hh1Ctl' - eta.*ca1ICtl(:,t);
        dc1Ctl = gm.*cc1Ctl' - om.*ca1CCtl(:,t);
    end
    
    ca3P(:,t+1) = ca3P(:,t) + da3;
    ca3I(:,t+1) = ca3I(:,t) + dh3;
    ca1PRmp(:,t+1) = ca1PRmp(:,t) + da1Rmp;
    ca1IRmp(:,t+1) = ca1IRmp(:,t) + dh1Rmp;
    ca1PSqr(:,t+1) = ca1PSqr(:,t) + da1Sqr;
    ca1ISqr(:,t+1) = ca1ISqr(:,t) + dh1Sqr;
    ca1PCtl(:,t+1) = ca1PCtl(:,t) + da1Ctl;
    ca1ICtl(:,t+1) = ca1ICtl(:,t) + dh1Ctl;
    
    if simTypeFlag == 2
        %Adaptation update
        ca3C(:,t+1)    = ca3C(:,t)    + dc3;
        ca1CRmp(:,t+1) = ca1CRmp(:,t) + dc1Rmp;
        ca1CSqr(:,t+1) = ca1CSqr(:,t) + dc1Sqr;
        ca1CCtl(:,t+1) = ca1CCtl(:,t) + dc1Ctl;
    end

end

actCell  = {ca3P,ca1PRmp,ca1PSqr,ca1PCtl};
hactCell = {ca3I,ca1IRmp,ca1ISqr,ca1ICtl};
inCell   = {ca3A,ca1ARmp,ca1ASqr,ca1ACtl};

%% Plot options

set(0,'DefaultLineLineWidth',2)

% Calculate Peaks
pks3  = zeros(N,1); locs3  = zeros(N,1);
pks1R = zeros(N,1); locs1R = zeros(N,1);
pks1S = zeros(N,1); locs1S = zeros(N,1);
pks1C = zeros(N,1); locs1C = zeros(N,1);

for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(ca3P(i,:));
    if ~isempty(locsTmp); pks3(i) = pksTmp(1); locs3(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(ca1PRmp(i,:));
    if ~isempty(locsTmp); pks1R(i) = pksTmp(1); locs1R(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(ca1PSqr(i,:));
    if ~isempty(locsTmp); pks1S(i) = pksTmp(1); locs1S(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(ca1PCtl(i,:));
    if ~isempty(locsTmp); pks1C(i) = pksTmp(1); locs1C(i) = locsTmp(1); end
end
pks1R = pks1R(pks1R>0); locs1R = locs1R(locs1R>0);
pks1S = pks1S(pks1S>0); locs1S = locs1S(locs1S>0);
locsCell = {locs3,locs1R,locs1S,locs1C};
pksCell = {pks3,pks1R,pks1S,pks1C};

pStruct.rampTypeFlag = rampTypeFlag;
pStruct.cueN = cueN;
pStruct.T = T;
pStruct.stimDelay = stimDelay;

% Threshold method
actThresh = 10;
tttCell = {};

for i = 1:numel(pksCell)
    pksTmp = pksCell{i};        %All the peaks for that sim (e.g. ramp)
    actTmp = actCell{i};        %Pyramidal activation for that sim e.g. NxT
    ttt = [];
    for j = 1:numel(pksTmp)     %For each peak in the sim
        indtmp = find(actTmp(j,:)>actThresh,1);
        if ~isempty(indtmp)
            ttt = [ttt,indtmp];
        end
    end
    tttCell(i) = {ttt};
end
if pStruct.rampTypeFlag ~= 'ctl'
    tttCell(2) = {tttCell{2}(tttCell{2} > stimDelay)};  %Only use the peaks after opto stim onset
    tttCell(3) = {tttCell{3}(tttCell{3} > stimDelay)};  %Only use the peaks after opto stim onset
end

ca1RmpPlot = plotCA1RipEx(ca3P,ca1PRmp,ca1ARmp,tttCell{2},pStruct);
set(ca1RmpPlot, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.6, 0.5, 0.38]);
% legend('CA3 Pyr','CA1 Ramp Pyr')
ca1SqrPlot = plotCA1RipEx(ca3P,ca1PSqr,ca1ASqr,tttCell{3},pStruct);
set(ca1SqrPlot, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.3, 0.5, 0.38]);
legend('CA3 Pyr','CA1 Square Pyr')
ca1CtlPlot = plotCA1RipEx(ca3P,ca1PCtl,ca1ACtl,tttCell{4},pStruct);
set(ca1CtlPlot, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.38]);
legend('CA3 Pyr','CA1 Control Pyr')

% % CA3 Pyr-CA3 Pyr Weight figure for manuscript
% figure(); baa=gca();
% imagesc(flipud(W'),[0,max(max(W))]); 
% title('CA3-CA3 Weight (W)');    
% xlabel('Pre-Synaptic CA3'); ylabel('Post-Synaptic CA3')
% baa.YTick = 1:7:15; baa.XTick = 1:7:15; 
% yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
% colorbar; axis square
% set(gca,'FontSize',20,'fontname','times')
% 
% CA3 Pyr-CA1 Pyr Weight figure for manuscript
figure();
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.2, 0.5, 0.8]);
% subplot(2,2,1); 
bad=gca();
imagesc(flipud(WZ'),[0,max(max(WZ))]); 
title('CA3 Pyr-CA1 Pyr Weight (WZ)');    
xlabel('Pre-Synaptic CA3'); ylabel('Post-Synaptic CA1')
bad.YTick = 1:7:15; bad.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',18,'fontname','times')

figure();
% subplot(2,2,2); 
bae=gca();
imagesc(flipud(QZ'),[0,max(max(QZ))]); 
title('CA1 IN-CA1 Pyr Weight (QZ)');    
xlabel('Pre-Synaptic CA1 IN'); ylabel('Post-Synaptic CA1 Pyr')
bae.YTick = 1:7:15; bae.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',18,'fontname','times')

figure();
% subplot(2,2,3); 
baf=gca();
imagesc(flipud(WQ'),[0,max(max(WQ))]); 
title('CA3 Pyr-CA1 IN Weight (WQ)');    
xlabel('Pre-Synaptic CA3 Pyr'); ylabel('Post-Synaptic CA1 IN')
baf.YTick = 1:7:15; baf.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',18,'fontname','times')

figure();
% subplot(2,2,4); 
baf=gca();
imagesc(flipud(ZQ'),[0,max(max(ZQ))]); 
title('CA1 Pyr-CA1 IN Weight (ZQ)');    
xlabel('Pre-Synaptic CA1 Pyr'); ylabel('Post-Synaptic CA1 IN')
baf.YTick = 1:7:15; baf.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',18,'fontname','times')

figure();
bag = gca; 
imagesc(flipud(ZZ'),[0,max(max(ZZ))]); 
title('CA1 Pyr-CA1 Pyr Weight (ZZ)');    
xlabel('Pre-Synaptic CA1 Pyr'); ylabel('Post-Synaptic CA1 Pyr')
bag.YTick = 1:7:15; bag.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',18,'fontname','times')

% CA3 IN-Pyr Weight figure for manuscript
figure(); bab=gca();
imagesc(flipud(H'),[0,max(max(H))]); 
% title('Pyr-Pyr Weight (W)');    
xlabel('Pre-Synaptic IN'); ylabel('Post-Synaptic CA3')
bab.YTick = 1:7:15; bab.XTick = 1:7:15; 
yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
colorbar; axis square
set(gca,'FontSize',20,'fontname','times')

% % Plot ANoise example
% figure;
% plot(ANoise(1,:),'k')
% % ylabel('Activity'); xlabel('Time')
% xlim([0 250]); ylim([-.15 .15])
% set(gca,'FontSize',20,'fontname','times')

%% Saves
if saveFlag == 1
disp('saving vars and figs')

sdel = num2str(stimDelay);

if rampTypeFlag == 'ctl'
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'fig')
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'png')
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'svg')
%     saveas(ca1SqrPlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'fig')
%     saveas(ca1SqrPlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'png')
%     saveas(ca1SqrPlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'svg')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'fig')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'png')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 1
    saveas(ca1RmpPlot,[saveDir,'FR_IMA_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'FR_IMA_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'FR_IMA_activity'],'svg')
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'fig')
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'png')
    saveas(ca1CtlPlot,[saveDir,'ctl_activity'],'svg')
    saveas(ca1SqrPlot,[saveDir,'sqr_activity'],'fig')
    saveas(ca1SqrPlot,[saveDir,'sqr_activity'],'png')
    saveas(ca1SqrPlot,[saveDir,'sqr_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 2
    saveas(ca1RmpPlot,[saveDir,'DR_IMA_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'DR_IMA_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'DR_IMA_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 5
    saveas(ca1RmpPlot,[saveDir,'BR_IMA_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'BR_IMA_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'BR_IMA_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 3
    saveas(ca1RmpPlot,[saveDir,'FR_IP_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'FR_IP_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'FR_IP_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 4
    saveas(ca1RmpPlot,[saveDir,'DR_IP_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'DR_IP_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'DR_IP_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 7
    saveas(ca1RmpPlot,[saveDir,'BR_IP_activity'],'fig')
    saveas(ca1RmpPlot,[saveDir,'BR_IP_activity'],'png')
    saveas(ca1RmpPlot,[saveDir,'BR_IP_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'svg')
end
disp('Done!')
end
