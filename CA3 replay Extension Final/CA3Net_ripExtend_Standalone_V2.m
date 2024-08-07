%%% Sample code for Replay Extension Model
%LKW 7/29/21
%Relies on plotRipExtendSingle.m for plots
%Calculates activity of CA3 pyr circuit undergoing a cued ripple/replay
%event with or without a secondary optogenetic pulse applied mid-replay. 
%Incorporates simple adaptation current using simulated intracellular Calcium dynamics
%Plots rate-based circuit activity. 
%%%

clearvars
close all
% rng(1)                          %For reproducibility of stochastic simulations

rampTypeFlag        = 12;        %'ctl' = no 2nd pulse; 1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 5 = BR IMA; 6 = Square Pulsatile IMA; 7 = BR IP; 8 = Square Pulsatile IP; 9 = Sinusoidal IMA; 10; Sinusoidal IP
simTypeFlag         = 2;        %1 = Linear; 2 = Linear with Adaptation
noiseFlag           = 0;        %0 = no noise; 1 = Voltage noise; 2 = Opsin Noise; 3 = Light Scattering noise; 4 = All noise
disp(['Ramp Type ', num2str(rampTypeFlag),'; Sim Type ', num2str(simTypeFlag), '; Noise type ', num2str(noiseFlag)]);
saveFlag            = 0;
saveDir = 'Your\Dir\Here\';

cueN                = 1;            %Cue node
N                   = 15;           %Nodes per region
Ww                  = 0.032;        %Weight strength pyr to pyr
Hh                  = 0.035;        %Weight strength IN to Pyr
Wh                  = 0.05;         %Weight strength pyr to IN
HAuto               = 0.003;        %Inhibition self feedback
tha                 = 4*ones(1,N);
thh                 = 4*ones(1,N);
thc                 = 4*ones(1,N);
eta                 = 0.01;         %Decay constant
T                   = 1500;         %Time steps, must be even
Iexcit1             = 1;            %Try 0.5 with Iexcit2 at 0.05 for ramp < square ripple
Iexcit2             = 0.09;         %Try as low as 0.09 or 0.05 opto pulse for good ripple spacing
inDur1              = 20;           %Duration for kicking off a ripple
inDur2              = 100;          %Stim duration
rampPerc            = 0.5;          %Percentage of ramp
rampLen             = round(inDur2*rampPerc);  %Duration ramp length e.g. 1/5 or 1/2
onsetDelay          = 50;           %Wait time to ripple start from sim start
stimDelay           = 150 + onsetDelay + inDur1;          %Wait time to opto pulse from simulation start
noiseAmp            = 0.5;          %Amplitude of Voltage noise. Try 0.1 or 0.5
noiseMu             = 1;            %Mean of ChR2 noise distro
noiseSigma          = 0.05;         %Variance of ChR2 noise distro

% Ionic Currents and related parameters 
Ek  = -10;              %Reversal Potential of Potassium
Ena = 70;               %Sodium
Ecl = 0;                %Chlorine
mu  = 0.01;             %Ca-dependent K-current
gm  = 0.001;            %Gamma; voltage-dependent Ca-currents
om  = 0.001;            %Omega; constant for diffusion of intracellular Ca

% Weight matrix
wtBias  = linspace(0.004,-0.004,N);
% wtBias  = zeros(1,N);
W       = zeros(N,N);    %Init wt mat
% W = rand(N,N).*Ww;          %All rand wt mat option
AH      = zeros(N,N);    %Wt mat of a to h (feedforward activation of IN)
H       = zeros(N,N);    %Weight of h to a (feedback inhibition)

for i = 1:N         %Build weight mats
    W(i,i)  = Ww+wtBias(i);   %Autorecurrency
    H(i,i)  = Hh;   %Direct feedback IN to Pyr
    AH(i,i) = Wh;   %Direct excitation pyr to IN
    if i <= N - 1   %Forward 1
%         W(i,i+1)  = (Ww+wtBias(i))/1.75;   %Pyr 2 Pyr
        W(i,i+1)  = (Ww+wtBias(i))/2;
%         H(i,i+1)  = Hh/2;   %IN  2 Pyr
%         AH(i,i+1) = Wh/2;   %Pyr 2 IN
    end
    if i <= N - 2  %Forward 2
%         W(i,i+2)  = (Ww+wtBias(i))/3.6;   %Pyr 2 Pyr
        W(i,i+2)  = (Ww+wtBias(i))/4;
%         H(i,i+2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i+2) = Wh/4;   %Pyr 2 IN
    end
    if i > 1    %Back 1
%         W(i,i-1)  = (Ww+wtBias(i))/1.75;   %Pyr 2 Pyr
%         W(i,i-1)  = Ww/2;   %Pyr 2 Pyr
%         H(i,i-1)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-1) = Wh/2;   %Pyr 2 IN
    end
    if i > 2    %Back 2   
%         W(i,i-2)  = (Ww+wtBias(i))/3.6;   %Pyr 2 Pyr
%         W(i,i-2)  = Ww/4;   %Pyr 2 Pyr
%         H(i,i-2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-2) = Wh/4;   %Pyr 2 IN
    end
end

%Ramp vectors
ARamp = zeros(N,T);   %Initialize input vector
aRamp = zeros(N,T);   %initialize a activity vector
hRamp = zeros(N,T);   %Initialize h activity vector
cRamp = zeros(N,T);   %Initialize Ca activity vector
%Square vectors
ASquare = zeros(N,T);
aSquare = zeros(N,T);
hSquare = zeros(N,T);
cSquare = zeros(N,T);
%Control vectors
AControl = zeros(N,T);
aControl = zeros(N,T);
hControl = zeros(N,T);
cControl = zeros(N,T);

%Set up waveform Afferent Inputs
ARamp(cueN,onsetDelay+1:onsetDelay+inDur1)   = Iexcit1;  %Start a ripple ramp
ASquare(cueN,onsetDelay+1:onsetDelay+inDur1) = Iexcit1;  %Start a ripple square
AControl(cueN,onsetDelay+1:onsetDelay+inDur1)= Iexcit1;  %Start a ripple control

ASquare(:,stimDelay+1:stimDelay+inDur2) = Iexcit2; % rand(N,inDur);  %Square pulse / random

if rampTypeFlag == 1        %FR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ARamp(:,stimDelay+1:stimDelay+rampLen)= repmat(linspace(0,Iexcit2,rampLen),[N,1]);   %Front Ramp
elseif rampTypeFlag == 2    %DR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2; 
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,Iexcit2,rampLen),[N,1]);  %Front Ramp
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
elseif rampTypeFlag == 5    %BR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
end
squarea = Iexcit2*inDur2;     %Area under the square pulse
if rampTypeFlag == 3        %FR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
elseif rampTypeFlag == 4    %DR IP
    IRamp = squarea/(inDur2 - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
elseif rampTypeFlag == 6    % Repeated Square Waves
    period = 20;           % Period of the square waves
    dutyCycle = 0.5;       % Duty cycle of the square waves
    pulseDuration = round(period * dutyCycle);
    numPulses = floor(inDur2 / period);
    for k = 0:numPulses-1
        startIdx = stimDelay + k * period + 1;
        endIdx = startIdx + pulseDuration - 1;
        ARamp(:,startIdx:endIdx) = Iexcit2;
    end
elseif rampTypeFlag == 7    %BR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
elseif rampTypeFlag == 8    % Iso-Power Repeated Square Waves
    period = 20;           % Period of the square waves
    dutyCycle = 0.5;       % Duty cycle of the square waves
    pulseDuration = round(period * dutyCycle);
    numPulses = floor(inDur2 / period);
    
    % Calculate the new pulse amplitude to maintain the same total area
    IRamp = squarea / (numPulses * pulseDuration);
    
    for k = 0:numPulses-1
        startIdx = stimDelay + k * period + 1;
        endIdx = startIdx + pulseDuration - 1;
        ARamp(:,startIdx:endIdx) = IRamp;
    end
elseif rampTypeFlag == 9    % Sinusoidal Waves with Same Max Amplitude as Square
    period = 20;           % Period of the sinusoidal waves
    numPulses = floor(inDur2 / period);
    
    % Generate time vector for one period
    t = linspace(0, period, period + 1);
    
    % Calculate the sinusoidal wave with half max amplitude of Iexcit2
    sinWave = (Iexcit2 / 2) * (sin(2 * pi * t / period - pi/2) + 1); % Sinusoid shifted up
    
    for k = 0:numPulses-1
        startIdx = stimDelay + k * period;
        endIdx = startIdx + period;
        ARamp(:,startIdx:endIdx) = repmat(sinWave, N, 1);
    end
elseif rampTypeFlag == 10    % Sinusoidal Waves with Same Area Under Curve as Square
    period = 20;           % Period of the sinusoidal waves
    numPulses = floor(inDur2 / period);
    
    % Generate time vector for one period
    t = linspace(0, period, period + 1);
    
    % Calculate the area of one sinusoidal wave period (shifted up)
    sinArea = trapz(t, sin(2 * pi * t / period - pi/2) + 1); % Shift the wave up so minimum is 0
    
    % Calculate the new amplitude to ensure the same total area
    IRamp = squarea / (numPulses * sinArea);
    
    % Generate the sinusoidal wave with the adjusted amplitude
    sinWave = IRamp * (sin(2 * pi * t / period - pi/2) + 1); % Adjust amplitude and shift up
    
    for k = 0:numPulses-1
        startIdx = stimDelay + k * period;
        endIdx = startIdx + period;
        ARamp(:,startIdx:endIdx) = repmat(sinWave, N, 1);
    end
elseif rampTypeFlag == 11 % Poisson Spike Train
    lambda = 5;
    pulseDuration = 5;
    refractoryPeriod = 2;

    n = 0;
    while n < inDur2
        n = n + max(poissrnd(lambda), refractoryPeriod);
        ARamp(:,stimDelay+n:stimDelay+n+pulseDuration)  = Iexcit2;
        n = n + pulseDuration
    end
elseif rampTypeFlag == 12    % Iso-Power Poisson Spike Train
    lambda = 5;    % Center frequency of the Poisson distribution (5-42 Hz)
    pulseDuration = 5;       % Pulse width in milliseconds (2-10 ms)
    refractoryPeriod = 2;
    % Initialize variables
    n = 0;                   % Start time
    pulseTimes = [];

    % Generate Poisson spike times
     while n < inDur2
        n = n + max(poissrnd(lambda), refractoryPeriod);
        pulseTimes = [pulseTimes, n];  % Store the pulse start time
        n = n + pulseDuration;
    end

    % Calculate the required pulse amplitude to match the total area
    numPulses = length(pulseTimes);
    squareArea = Iexcit2 * inDur2;  % Total area under the square pulse
    pulseArea = pulseDuration * numPulses;  % Total duration of all pulses
    IRamp = squareArea / pulseArea;  % Adjusted amplitude to match total area

    % Apply pulses to ARamp
    for k = 1:numPulses
        startIdx = stimDelay + pulseTimes(k);  % +1 to convert to 1-based index
        endIdx = startIdx + pulseDuration;
        if endIdx <= T
            ARamp(:, startIdx:endIdx) = IRamp;  % Apply adjusted amplitude
        end
    end
end



%% Build noise

if noiseFlag == 0
    ANoise = zeros(N,T);
    tissNoise = ones(N,1);
elseif noiseFlag == 1   %V Noise only
%     rng(kern);                  %Noise kernel on/off
    ChR2Noise = ones(N,1);        %For testing V Noise alone
    distNoise = ones(N,1);        %For testing V Noise alone
    ANoise = 2*noiseAmp.*(rand(N,T) - 0.5);
    
    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect

elseif noiseFlag == 2   %ChR2 Noise only
    ANoise = zeros(N,T);
    distNoise = ones(N,1);        %For testing ChR2Noise alone

    ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression

    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect
    
elseif noiseFlag == 3   %Distance Noise only
    ANoise = zeros(N,T);
    ChR2Noise = ones(N,1);
    
%     Light Scattering Noise from distance in 3D space
%     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
    load ScatteringFitIrr8mW.mat    %fspline calculated from 8mW irradiance (mW/mm^2)
%     load ScatteringFitIrr10mW.mat    %fspline calculated from 10mW irradiance (mW/mm^2)
    sTip = [0 0 0]; %set laser source at origin
    unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
    unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
    distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
    distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
    distNoise = distNoise/5;        %Get normalized % activation with distance dropoff

    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect

elseif noiseFlag == 4   %Combined Noise
    ANoise = 2*noiseAmp.*(rand(N,T) - 0.5);

    ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression

%     Light Scattering Noise from distance in 3D space
%     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
%     load ScatteringFitIrr8mW.mat    %fspline calculated from 8mW irradiance (mW/mm^2)
    load ScatteringFitIrr10mW.mat    %fspline calculated from 10mW irradiance (mW/mm^2)
    sTip = [0 0 0]; %set laser source at origin
    unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
    unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
    distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
    distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
    distNoise = distNoise/5;        %Get normalized % activation with distance dropoff
        
    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect
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
        if (aRamp(j,t)   - tha(j)) > 0; aaRamp(j)   = aRamp(j,t)   - tha(j); else aaRamp(j)   = 0; end %Threshold activity over 'tha' to 0
        if (aSquare(j,t) - tha(j)) > 0; aaSquare(j) = aSquare(j,t) - tha(j); else aaSquare(j) = 0; end
        if (aControl(j,t)- tha(j)) > 0; aaControl(j)= aControl(j,t)- tha(j); else aaControl(j)= 0; end
        if (hRamp(j,t)   - thh(j)) > 0; hhRamp(j)   = hRamp(j,t)   - thh(j); else hhRamp(j)   = 0; end %Threshold activity over 'thh' to 0
        if (hSquare(j,t) - thh(j)) > 0; hhSquare(j) = hSquare(j,t) - thh(j); else hhSquare(j) = 0; end
        if (hControl(j,t)- thh(j)) > 0; hhControl(j)= hControl(j,t)- thh(j); else hhControl(j)= 0; end
        if (aRamp(j,t)   - thc(j)) > 0; ccRamp(j)   = aRamp(j,t)   - thc(j); else ccRamp(j)   = 0; end %Threshold activity under 'thc' to 0
        if (aSquare(j,t) - thc(j)) > 0; ccSquare(j) = aSquare(j,t) - thc(j); else ccSquare(j) = 0; end
        if (aControl(j,t)- thc(j)) > 0; ccControl(j)= aControl(j,t)- thc(j); else ccControl(j)= 0; end
    end
    
    if simTypeFlag == 1
        %Standard Linear without Adapatation
        daRamp = tissNoise.*ARamp(:,t) + ANoise(:,t) + (aaRamp*W)' - (hhRamp*H)' - eta.*aRamp(:,t);
        dhRamp = (aaRamp*AH)' - HAuto.*hhRamp' - eta.*hRamp(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        daSquare = tissNoise.*ASquare(:,t) + ANoise(:,t) + (aaSquare*W)' - (hhSquare*H)' - eta.*aSquare(:,t);
        dhSquare = (aaSquare*AH)' - HAuto.*hhSquare' - eta.*hSquare(:,t);
        daControl = tissNoise.*AControl(:,t) + ANoise(:,t) + (aaControl*W)' - (hhControl*H)' - eta.*aControl(:,t);
        dhControl = (aaControl*AH)' - HAuto.*hhControl' - eta.*hControl(:,t);
    elseif simTypeFlag == 2
        %Linear variant with Adaptation
        daRamp = tissNoise.*ARamp(:,t) + ANoise(:,t) + (aaRamp*W)' - (hhRamp*H)' - eta.*aRamp(:,t) + mu.*cRamp(:,t).*(Ek - aRamp(:,t));
        dhRamp = (aaRamp*AH)' - HAuto.*hhRamp' - eta.*hRamp(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        dcRamp = gm.*ccRamp' - om.*cRamp(:,t);
        daSquare = tissNoise.*ASquare(:,t) + ANoise(:,t) + (aaSquare*W)' - (hhSquare*H)' - eta.*aSquare(:,t) + mu.*cSquare(:,t).*(Ek - aSquare(:,t));
        dhSquare = (aaSquare*AH)' - HAuto.*hhSquare' - eta.*hSquare(:,t);
        dcSquare = gm.*ccSquare' - om.*cSquare(:,t);
        daControl = tissNoise.*AControl(:,t) + ANoise(:,t) + (aaControl*W)' - (hhControl*H)' - eta.*aControl(:,t) + mu.*cControl(:,t).*(Ek - aControl(:,t));
        dhControl = (aaControl*AH)' - HAuto.*hhControl' - eta.*hControl(:,t);
        dcControl = gm.*ccControl' - om.*cControl(:,t);
    end
    
    aRamp(:,t+1) = aRamp(:,t) + daRamp; %Update next time point of a
    hRamp(:,t+1) = hRamp(:,t) + dhRamp; %Update next h time point
    aSquare(:,t+1) = aSquare(:,t) + daSquare;
    hSquare(:,t+1) = hSquare(:,t) + dhSquare;
    aControl(:,t+1) = aControl(:,t) + daControl;
    hControl(:,t+1) = hControl(:,t) + dhControl;
    
    if simTypeFlag == 2
        %Adaptation update
        cRamp(:,t+1) = cRamp(:,t) + dcRamp;
        cSquare(:,t+1) = cSquare(:,t) + dcSquare;
        cControl(:,t+1) = cControl(:,t) + dcControl;
    end

end

actCell = {aRamp,aSquare,aControl};
hactCell = {hRamp,hSquare,hControl};
inCell = {ARamp,ASquare,AControl};

%% Plot options

set(0,'DefaultLineLineWidth',2)

% Calculate Peaks
pksR = zeros(N,1); locsR = zeros(N,1);
pksS = zeros(N,1); locsS = zeros(N,1);
pksC = zeros(N,1); locsC = zeros(N,1);

for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(aRamp(i,:));
    if ~isempty(locsTmp); pksR(i) = pksTmp(1); locsR(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(aSquare(i,:));
    if ~isempty(locsTmp); pksS(i) = pksTmp(1); locsS(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(aControl(i,:));
    if ~isempty(locsTmp); pksC(i) = pksTmp(1); locsC(i) = locsTmp(1); end
end
pksC = pksC(pksC>0); locsC = locsC(locsC>0);
locsCell = {locsR,locsS,locsC};
pksCell = {pksR,pksS,pksC};

pStruct.rampTypeFlag = rampTypeFlag;
pStruct.cueN = cueN;
pStruct.T = T;
pStruct.stimDelay = stimDelay;

% Threshold method
actThresh = 10;
tttCell = {};

for i = 1:3
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
    tttCell(1) = {tttCell{1}(tttCell{1} > stimDelay)};  %Only use the peaks after opto stim onset
end
rampPlot = plotRipExtendSingle_V2(aRamp,ARamp,tttCell{1},pStruct);
rampINPlot = plotRipExtend_IN_V2(aRamp,hRamp,ARamp,pStruct);
if rampTypeFlag == 'ctl'
    pStruct.rampTypeFlag = 'sqr';
    squarePlot = plotRipExtendSingle_V2(aSquare,ASquare,tttCell{2},pStruct);
    squareINPlot = plotRipExtend_IN_V2(aSquare,hSquare,ASquare,pStruct);
end

% % Pyr-Pyr Weight figure for manuscript
% figure(); baa=gca();
% imagesc(flipud(W'),[0,max(max(W))]); 
% % title('Pyr-Pyr Weight (W)');    
% xlabel('Pre-Synaptic CA3'); ylabel('Post-Synaptic CA3')
% baa.YTick = 1:7:15; baa.XTick = 1:7:15; 
% yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
% colorbar; axis square
% set(gca,'FontSize',20,'fontname','times')

% % IN-Pyr Weight figure for manuscript
% figure(); bab=gca();
% imagesc(flipud(H'),[0,max(max(H))]); 
% % title('Pyr-Pyr Weight (W)');    
% xlabel('Pre-Synaptic IN'); ylabel('Post-Synaptic CA3')
% bab.YTick = 1:7:15; bab.XTick = 1:7:15; 
% yticklabels(linspace(15,1,3)); xticklabels(linspace(1,15,3));
% colorbar; axis square
% set(gca,'FontSize',20,'fontname','times')

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
    saveas(rampPlot,[saveDir,'ctl_activity'],'fig')
    saveas(rampPlot,[saveDir,'ctl_activity'],'png')
    saveas(rampPlot,[saveDir,'ctl_activity'],'svg')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'ctl_IN_activity'],'svg')
    saveas(squarePlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'fig')
    saveas(squarePlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'png')
    saveas(squarePlot,[saveDir,'sqr_activity','_stimDelay_',sdel],'svg')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'fig')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'png')
%     saveas(squareINPlot,[saveDir,'sqr_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 1
    saveas(rampPlot,[saveDir,'FR_IMA_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'FR_IMA_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'FR_IMA_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'FR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 2
    saveas(rampPlot,[saveDir,'DR_IMA_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'DR_IMA_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'DR_IMA_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'DR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 5
    saveas(rampPlot,[saveDir,'BR_IMA_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'BR_IMA_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'BR_IMA_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IMA_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 3
    saveas(rampPlot,[saveDir,'FR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'FR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'FR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'FR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 4
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 6
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'DR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'DR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 7
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 8
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 9
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'svg')
elseif pStruct.rampTypeFlag == 10
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'fig')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'png')
    saveas(rampPlot,[saveDir,'BR_IP_activity','_stimDelay_',sdel],'svg')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'fig')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'png')
%     saveas(rampINPlot,[saveDir,'BR_IP_IN_activity'],'svg')
end
disp('Saving done')
end
