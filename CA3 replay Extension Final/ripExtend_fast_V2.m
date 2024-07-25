function [actCell,hactCell,inCell] = ripExtend_fast_V2(pStruct)
%LKW 8/2/21
%Inputs: 
%pStruct = A struct containing fields set up in the CA3Net_ripExtend_Search_V2.m function or another script
%Outputs: 
% actCell = cell of matrices containing the estimated activity of Pyrs for each manipulation (ramp, square, control = no opto)
% hactCell = cell of matrices containing the estimated activity of INs for each manipulation
% inCell = cell of matrices containing the input delivered for each manipulation

%Unpack pStruct
cueN    = pStruct.cueN;         %Node for starting the ripple
N       = pStruct.N;            %Number of Nodes
T       = pStruct.T;            %Lengh of Time Series
tha     = pStruct.tha;          %Vector of Pyr thresholds
thh     = pStruct.thh;          %Vector of IN thresholds
eta     = pStruct.eta;          %Leak current of all nodes
HAuto   = pStruct.HAuto;        %IN self inhibition
thc     = pStruct.thc; 
mu      = pStruct.mu; 
gm      = pStruct.gm; 
om      = pStruct.om;
Ek      = -10;                  %Reversal potentials of Potassium
Iexcit1 = pStruct.Iexcit1;      %Strength of cue impulse
Iexcit2 = pStruct.Iexcit2;      %Strength of opto impulse
inDur1  = pStruct.inDur1;       %Length of cue impulse
inDur2  = round(pStruct.inDur2);%Length of opto impulse

%Set up Weight matrices
W = pStruct.W; H = pStruct.H; AH = pStruct.AH;

%Make new set of vectors to store pyramidals, INs, and inputs
%Ramp vectors
ARamp = zeros(N,T);   %Initialize input vector
aRamp = zeros(N,T);   %initialize a activity vector
hRamp = zeros(N,T);   %Initialize h activity vector
cRamp = zeros(N,T);   % c activity vector
%Control vectors
AControl = zeros(N,T);
aControl = zeros(N,T);
hControl = zeros(N,T);
cControl = zeros(N,T);

%Initialize inputs for afferent input A
onsetDelay         = pStruct.onsetDelay;        %Wait time to ripple start
stimDelay          = round(pStruct.stimDelay);  %Wait time to opto pulse
rampLen            = pStruct.rampLen;           %Duration ramp length

ARamp(cueN,onsetDelay+1:onsetDelay+inDur1)   = Iexcit1;  %Start a ripple ramp
AControl(cueN,onsetDelay+1:onsetDelay+inDur1)= Iexcit1;  %Start a ripple control

% Setup afferent pulses
if pStruct.rampTypeFlag == 1        %FR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ARamp(:,stimDelay+1:stimDelay+rampLen)= repmat(linspace(0,Iexcit2,rampLen),[N,1]);   %Front Ramp
elseif pStruct.rampTypeFlag == 2    %DR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2; 
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,Iexcit2,rampLen),[N,1]);  %Front Ramp
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
elseif pStruct.rampTypeFlag == 5    %BR IMA
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
elseif pStruct.rampTypeFlag == 6    % Evenly spaced pulses
    dutycycle = 0.4;
    numPulse = 20;
    period = inDur2 / numPulse;
    pulseLength = period * dutycycle;
    rampLength = floor(pStruct.tmpRamp * pulseLength);
    for i = 1:numPulse
        % Square
        % ARamp(:,floor(stimDelay+(i-1)*period):floor(stimDelay+(i-1)*period+pulseLength)) = Iexcit2;

        % FR IMA

        pulseStart = floor(stimDelay + period * (i - 1));
        pulseEnd = floor(pulseStart + pulseLength);
        rampEnd = floor(pulseStart + rampLength);

        ARamp(:,pulseStart:pulseEnd)  = Iexcit2;
        % disp(linspace(0,Iexcit2,rampLen));
        ARamp(:,pulseStart:rampEnd) = repmat(linspace(0,Iexcit2,rampEnd - pulseStart + 1),[N,1]);
    end
elseif pStruct.rampTypeFlag == 8 % Sinusoidal
    % frequency = 30; % hz
    frequency = pStruct.tmpFrequency;
    period = 1000 / frequency;

    ARamp(:,stimDelay+1:stimDelay+1+inDur2) = (sin(repmat(linspace(0,inDur2,inDur2 + 1)*2*pi/period, [N,1])) * Iexcit2 + Iexcit2) / 2;
elseif pStruct.rampTypeFlag == 9 % Poisson Spike Train
    centerFrequency = 10;
    pulseDuration = 5;
    refractoryPeriod = 1;

    n = 1;
    while n < inDur2
        n = max(n + poissrnd(centerFrequency), n+pulseDuration+refractoryPeriod);
        ARamp(:,stimDelay+n+1:stimDelay+1+n+pulseDuration)  = Iexcit2;
    end

end



squarea = Iexcit2*inDur2;     %Area under the square pulse
if pStruct.rampTypeFlag == 3        %FR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
elseif pStruct.rampTypeFlag == 4    %DR IP
    IRamp = squarea/(inDur2 - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
elseif pStruct.rampTypeFlag == 7    %BR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ARamp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
end

%Build noise
if pStruct.noiseFlag == 3 || pStruct.noiseFlag == 4
if pStruct.IrrVal == 8
    load ScatteringFitIrr8mW.mat    %fspline calculated from 8mW irradiance (mW/mm^2)
elseif pStruct.IrrVal == 10
    load ScatteringFitIrr10mW.mat    %fspline calculated from 10mW irradiance (mW/mm^2)
end
end

if pStruct.noiseFlag == 0
    ANoise = zeros(N,T);
    tissNoise = ones(N,1);
elseif pStruct.noiseFlag == 1   %V Noise only
    noiseAmp = pStruct.noiseAmp;
%     rng(pStruct.kern);                  %Noise kernel on/off
    ChR2Noise = ones(N,1);        %For testing V Noise alone
    distNoise = ones(N,1);        %For testing V Noise alone
    ANoise = 2*noiseAmp.*(rand(N,T) - 0.5);
    
    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect

elseif pStruct.noiseFlag == 2   %ChR2 Noise only
    ANoise = zeros(N,T);
    distNoise = ones(N,1);        %For testing ChR2Noise alone

    noiseMu = pStruct.noiseMu;
    noiseSigma = pStruct.noiseSigma;
    ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression

    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect
    
elseif pStruct.noiseFlag == 3   %Distance Noise only
    ANoise = zeros(N,T);
    ChR2Noise = ones(N,1);
    
%     Light Scattering Noise from distance in 3D space
%     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
    sTip = [0 0 0]; %set laser source at origin
    unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
    unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
    distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
    distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
    distNoise = distNoise/5;        %Get normalized % activation with distance dropoff

    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect

elseif pStruct.noiseFlag == 4   %Combined Noise
    noiseAmp = pStruct.noiseAmp;
    ANoise = 2*noiseAmp.*(rand(N,T) - 0.5);

    noiseMu = pStruct.noiseMu;
    noiseSigma = pStruct.noiseSigma;
    ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression

%     Light Scattering Noise from distance in 3D space
%     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
    sTip = [0 0 0]; %set laser source at origin
    unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
    unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
    distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
    distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
    distNoise = distNoise/5;        %Get normalized % activation with distance dropoff
        
    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect
end

% Run simulation
for t=1:T-1
    for j = 1:N
        if (aRamp(j,t)   - tha(j)) > 0; aaRamp(j)   = aRamp(j,t)   - tha(j); else aaRamp(j)   = 0; end %Threshold activity over 'tha' to 0
        if (aControl(j,t)- tha(j)) > 0; aaControl(j)= aControl(j,t)- tha(j); else aaControl(j)= 0; end
        if (hRamp(j,t)   - thh(j)) > 0; hhRamp(j)   = hRamp(j,t)   - thh(j); else hhRamp(j)   = 0; end %Threshold activity over 'thh' to 0
        if (hControl(j,t)- thh(j)) > 0; hhControl(j)= hControl(j,t)- thh(j); else hhControl(j)= 0; end
        if (aRamp(j,t)   - thc(j)) > 0; ccRamp(j)   = aRamp(j,t)   - thc(j); else ccRamp(j)   = 0; end %Threshold activity under 'thc' to 0
        if (aControl(j,t)- thc(j)) > 0; ccControl(j)= aControl(j,t)- thc(j); else ccControl(j)= 0; end
    end
        
    if pStruct.simTypeFlag == 1
        %Standard Linear without Adapatation
        daRamp = tissNoise.*ARamp(:,t) + ANoise(:,t) + (aaRamp*W)' - (hhRamp*H)' - eta.*aRamp(:,t);
        dhRamp = (aaRamp*AH)' - HAuto.*hhRamp' - eta.*hRamp(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        daControl = tissNoise.*AControl(:,t) + ANoise(:,t) + (aaControl*W)' - (hhControl*H)' - eta.*aControl(:,t);
        dhControl = (aaControl*AH)' - HAuto.*hhControl' - eta.*hControl(:,t);
    elseif pStruct.simTypeFlag == 2
        %Linear variant with Adaptation
        daRamp = tissNoise.*ARamp(:,t) + ANoise(:,t) + (aaRamp*W)' - (hhRamp*H)' - eta.*aRamp(:,t) + mu.*cRamp(:,t).*(Ek - aRamp(:,t));
        dhRamp = (aaRamp*AH)' - HAuto.*hhRamp' - eta.*hRamp(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        dcRamp = gm.*ccRamp' - om.*cRamp(:,t);
        daControl = tissNoise.*AControl(:,t) + ANoise(:,t) + (aaControl*W)' - (hhControl*H)' - eta.*aControl(:,t) + mu.*cControl(:,t).*(Ek - aControl(:,t));
        dhControl = (aaControl*AH)' - HAuto.*hhControl' - eta.*hControl(:,t);
        dcControl = gm.*ccControl' - om.*cControl(:,t);
    end
    
    aRamp(:,t+1) = aRamp(:,t) + daRamp; %Update next time point of a
    hRamp(:,t+1) = hRamp(:,t) + dhRamp; %Update next h time point
    aControl(:,t+1) = aControl(:,t) + daControl;
    hControl(:,t+1) = hControl(:,t) + dhControl;
    
    if pStruct.simTypeFlag == 2
        %Adaptation update
        cRamp(:,t+1) = cRamp(:,t) + dcRamp;
        cControl(:,t+1) = cControl(:,t) + dcControl;
    end

end

% Format outputs
actCell     = {aRamp,aControl};
hactCell    = {hRamp,hControl};
inCell      = {ARamp,AControl};

end