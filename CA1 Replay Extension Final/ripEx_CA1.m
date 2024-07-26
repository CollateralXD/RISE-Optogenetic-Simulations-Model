function [actCell,hactCell,inCell] = ripEx_CA1(pStruct)
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
W = pStruct.W; H = pStruct.H; AH = pStruct.AH; ZZ = pStruct.ZZ;
WZ = pStruct.WZ; WQ = pStruct.WQ; QZ = pStruct.QZ; ZQ = pStruct.ZQ;

%Make new set of vectors to store pyramidals, INs, and inputs
%CA3 vectors
ca3A = zeros(N,T);
ca3P = zeros(N,T);
ca3I = zeros(N,T);
ca3C = zeros(N,T);
%CA1 Ramp vectors
ca1ARmp = zeros(N,T);   %Initialize input vector
ca1PRmp = zeros(N,T);   %initialize a activity vector
ca1IRmp = zeros(N,T);   %Initialize h activity vector
ca1CRmp = zeros(N,T);   % c activity vector
%CA1 Control vectors
ca1ACtl = zeros(N,T);
ca1PCtl = zeros(N,T);
ca1ICtl = zeros(N,T);
ca1CCtl = zeros(N,T);

%Initialize inputs for afferent input A
onsetDelay         = pStruct.onsetDelay;        %Wait time to ripple start
stimDelay          = round(pStruct.stimDelay);  %Wait time to opto pulse
rampLen            = pStruct.rampLen;           %Duration ramp length

% ca1ARmp(cueN,onsetDelay+1:onsetDelay+inDur1)   = Iexcit1;  %Start a ripple ramp
% ca1ACtl(cueN,onsetDelay+1:onsetDelay+inDur1)= Iexcit1;  %Start a ripple control
ca3A(cueN,onsetDelay+1:onsetDelay+inDur1) = Iexcit1;    %Start CA3 ripple

% Setup afferent pulses
if pStruct.rampTypeFlag == 1        %FR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen)= repmat(linspace(0,Iexcit2,rampLen),[N,1]);   %Front Ramp
elseif pStruct.rampTypeFlag == 2    %DR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2; 
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,Iexcit2,rampLen),[N,1]);  %Front Ramp
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
elseif pStruct.rampTypeFlag == 5    %BR IMA
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = Iexcit2;
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(Iexcit2,0,rampLen),[N,1]);  %Rear ramp
end
squarea = Iexcit2*inDur2;     %Area under the square pulse
if pStruct.rampTypeFlag == 3        %FR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
elseif pStruct.rampTypeFlag == 4    %DR IP
    IRamp = squarea/(inDur2 - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+1:stimDelay+rampLen) = repmat(linspace(0,IRamp,rampLen),[N,1]);   %Front Ramp
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
elseif pStruct.rampTypeFlag == 7    %BR IP
    IRamp = squarea/(rampLen/2 + inDur2 - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ca1ARmp(:,stimDelay+1:stimDelay+inDur2)  = IRamp;
    ca1ARmp(:,stimDelay+inDur2-rampLen+1:stimDelay+inDur2) = repmat(linspace(IRamp,0,rampLen),[N,1]);   %Rear ramp
end

%Build noise
if pStruct.noiseFlag == 0
    ANoise = zeros(N,T);
    tissNoise = ones(N,1);
elseif pStruct.noiseFlag == 1
    noiseAmp = pStruct.noiseAmp;
    noiseMu = pStruct.noiseMu;
    noiseSigma = pStruct.noiseSigma;
%     rng(pStruct.kern);                  %Noise kernel on/off

    %Voltage Noise from LFP fluctuations
    ANoise = noiseAmp.*(rand(N,T) - 0.5);
%     ChR2Noise = ones(N,1);        %For testing ANoise alone
%     distNoise = ones(N,1);        %For testing ANoise alone

% %     Light Scattering Noise from distance in 3D space
% %     load ScatteringFitVars.mat      %fspline calculated from percent power dropoff with distance through tissue
% %     load ScatteringFitIrr8mW.mat    %fspline calculated from 8mW irradiance (mW/mm^2)
%     load ScatteringFitIrr10mW.mat    %fspline calculated from 10mW irradiance (mW/mm^2)
%     sTip = [0 0 0]; %set laser source at origin
%     unitLocs = [rand(1,N)-0.5;rand(1,N)-0.5;rand(1,N)/10+0.2]; % Arranged 3 x N where rows are X, Y, Z in mm
%     unitDist = sqrt((sTip(1) - unitLocs(1,:)).^2 + (sTip(2) - unitLocs(2,:)).^2 + (sTip(3) - unitLocs(3,:)).^2);
%     distNoise = fspline(unitDist);  %Apply distance to power dropoff transform
%     distNoise(distNoise > 5) = 5;   %Threshold irradiance > 5 to 5mW
%     distNoise = distNoise/5;        %Get normalized % activation with distance dropoff
%     ChR2Noise = ones(N,1);        %For testing distNoise alone
    
    %Protein Expression Noise
    ChR2Noise = noiseMu - abs(randn(N,1)*noiseSigma); %Gain factor applied to each afferent pulse reflecting heterogeneity of ChR2 expression
    distNoise = ones(N,1);        %For testing ChR2Noise alone
    
    tissNoise = ChR2Noise .* distNoise;   %Combined noise effect
end

% Run simulation
for t=1:T-1
    for j = 1:N
        %CA3 Region
        if (ca3P(j,t)    - tha(j)) > 0; aa3(j)    = ca3P(j,t)    - tha(j); else aa3(j)    = 0; end %Threshold activity over 'tha' to 0
        if (ca3I(j,t)    - thh(j)) > 0; hh3(j)    = ca3I(j,t)    - thh(j); else hh3(j)    = 0; end %Threshold activity over 'thh' to 0
        if (ca3P(j,t)    - thc(j)) > 0; cc3(j)    = ca3P(j,t)    - thc(j); else cc3(j)    = 0; end %Threshold activity under 'thc' to 0
        %CA1 Region
        if (ca1PRmp(j,t) - tha(j)) > 0; aa1Rmp(j) = ca1PRmp(j,t) - tha(j); else aa1Rmp(j) = 0; end %Threshold activity over 'tha' to 0
        if (ca1PCtl(j,t) - tha(j)) > 0; aa1Ctl(j) = ca1PCtl(j,t) - tha(j); else aa1Ctl(j) = 0; end
        if (ca1IRmp(j,t) - thh(j)) > 0; hh1Rmp(j) = ca1IRmp(j,t) - thh(j); else hh1Rmp(j) = 0; end %Threshold activity over 'thh' to 0
        if (ca1ICtl(j,t) - thh(j)) > 0; hh1Ctl(j) = ca1ICtl(j,t) - thh(j); else hh1Ctl(j) = 0; end
        if (ca1PRmp(j,t) - thc(j)) > 0; cc1Rmp(j) = ca1PRmp(j,t) - thc(j); else cc1Rmp(j) = 0; end %Threshold activity under 'thc' to 0
        if (ca1PCtl(j,t) - thc(j)) > 0; cc1Ctl(j) = ca1PCtl(j,t) - thc(j); else cc1Ctl(j) = 0; end
    end
    
    if pStruct.simTypeFlag == 1
        %Standard Linear without Adapatation
        %CA3
        da3    = tissNoise.*ca3A(:,t) + ANoise(:,t) + (aa3*W)' - (hh3*H)' - eta.*ca3P(:,t);
        dh3    = (aa3*AH)' - HAuto.*hh3' - eta.*ca3I(:,t);
        %CA1
        da1Rmp = tissNoise.*ca1ARmp(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Rmp*ZZ)' - (hh1Rmp*QZ)' - eta.*ca1PRmp(:,t);
        dh1Rmp = (aa3*WQ)' + (aa1Rmp*ZQ)' - HAuto.*hh1Rmp' - eta.*ca1IRmp(:,t);
        da1Ctl = tissNoise.*ca1ACtl(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Ctl*ZZ)' - (hh1Ctl*QZ)' - eta.*ca1PCtl(:,t);
        dh1Ctl = (aa3*WQ)' + (aa1Ctl*ZQ)' - HAuto.*hh1Ctl' - eta.*ca1ICtl(:,t);
    elseif pStruct.simTypeFlag == 2
        da3    = tissNoise.*ca3A(:,t) + ANoise(:,t) + (aa3*W)' - (hh3*H)' - eta.*ca3P(:,t) + mu.*ca3C(:,t).*(Ek - ca3P(:,t));
        dh3    = (aa3*AH)' - HAuto.*hh3' - eta.*ca3I(:,t);
        dc3    = gm.*cc3' - om.*ca3C(:,t);
        %Linear variant with Adaptation
        da1Rmp = tissNoise.*ca1ARmp(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Rmp*ZZ)' - (hh1Rmp*QZ)' - eta.*ca1PRmp(:,t) + mu.*ca1CRmp(:,t).*(Ek - ca1PRmp(:,t));
        dh1Rmp = (aa3*WQ)' + (aa1Rmp*ZQ)' - HAuto.*hh1Rmp' - eta.*ca1IRmp(:,t);
        dc1Rmp = gm.*cc1Rmp' - om.*ca1CRmp(:,t);
        da1Ctl = tissNoise.*ca1ACtl(:,t) + ANoise(:,t) + (aa3*WZ)' + (aa1Ctl*ZZ)' - (hh1Ctl*QZ)' - eta.*ca1PCtl(:,t) + mu.*ca1CCtl(:,t).*(Ek - ca1PCtl(:,t));
        dh1Ctl = (aa3*WQ)' + (aa1Ctl*ZQ)' - HAuto.*hh1Ctl' - eta.*ca1ICtl(:,t);
        dc1Ctl = gm.*cc1Ctl' - om.*ca1CCtl(:,t);
    end
    
    ca3P(:,t+1) = ca3P(:,t) + da3;
    ca3I(:,t+1) = ca3I(:,t) + dh3;
    ca1PRmp(:,t+1) = ca1PRmp(:,t) + da1Rmp;
    ca1IRmp(:,t+1) = ca1IRmp(:,t) + dh1Rmp;
    ca1PCtl(:,t+1) = ca1PCtl(:,t) + da1Ctl;
    ca1ICtl(:,t+1) = ca1ICtl(:,t) + dh1Ctl;
    
    if pStruct.simTypeFlag == 2
        %Adaptation update
        ca3C(:,t+1)    = ca3C(:,t)    + dc3;
        ca1CRmp(:,t+1) = ca1CRmp(:,t) + dc1Rmp;
        ca1CCtl(:,t+1) = ca1CCtl(:,t) + dc1Ctl;
    end

end

% Format outputs
actCell     = {ca3P,ca1PRmp,ca1PCtl};
hactCell    = {ca3I,ca1IRmp,ca1ICtl};
inCell      = {ca3A,ca1ARmp,ca1ACtl};

end