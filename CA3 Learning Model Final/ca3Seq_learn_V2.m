function [actCell,hactCell,inCell,wOut] = ca3Seq_learn_V2(pStruct)
%Inputs: 
%pStruct = A struct containing fields set up in the CA3Net_learn_Search.m function or another script
%Outputs: 
% actCell = cell of matrices containing the estimated activity of Pyrs for each manipulation (ramp, square, control = no opto)
% hactCell = cell of matrices containing the estimated activity of INs for each manipulation
% inCell = cell of matrices containing the input delivered for each manipulation
% wOut   = cell of final weight snap shot for this simulation
%LKW 8/17/21

% rng(1)
%Unpack pStruct
cueN    = pStruct.cueN;         %Node for starting the ripple
dt      = pStruct.dt;
N       = pStruct.N;            %Number of Nodes
T       = pStruct.TLearn;       %Number of learn dts
T2      = pStruct.TTest;        %Number of test dts
tha     = pStruct.tha;          %Vector of Pyr thresholds
thh     = pStruct.thh;          %Vector of IN thresholds
eta1    = pStruct.eta1;        %Leak current of all nodes
eta2    = pStruct.eta2;         %Just in case
HAuto   = pStruct.HAuto;        %IN self inhibition
achLearn= pStruct.achLearn;     %Ach state of autorecurrency
achTest = pStruct.achTest;
Ww      = pStruct.Ww;
e1      = pStruct.e1;
e3      = pStruct.e3;
thc     = pStruct.thc; 
mu      = pStruct.mu; 
gm      = pStruct.gm; 
om      = pStruct.om;
Ek      = -10;                  %Reversal potentials of Potassium

Iexcit1 = pStruct.Iexcit1;      %Strength of learn pulses
Iexcit2 = pStruct.Iexcit2;      %Strength of test pulses
learnDur= round(pStruct.learnDur);     %Length of cue pulse
testDur = pStruct.testDur;      %Length of opto impulse

%Set up Weight matrices
H = pStruct.H; AH = pStruct.AH;
if pStruct.wtNoiseFlag == 0
    W  = pStruct.W;
elseif pStruct.wtNoiseFlag == 1
    W = lognrnd(pStruct.noiseMu,pStruct.noiseSigma,N,N).*pStruct.noiseAmp;          %Lognormal rand wt mat
end
W_pre = W;

%Make new set of vectors to store pyramidals, INs, and inputs
%Learn vectors
AIn1    = zeros(N,T);   %Initialize input vector
aLearn  = zeros(N,T);   %initialize a activity vector
hLearn  = zeros(N,T);   %Initialize h activity vector
cLearn  = zeros(N,T);   % c activity vector
%Test vectors
AIn2    = zeros(N,T);
aTest   = zeros(N,T);
hTest   = zeros(N,T);
cTest   = zeros(N,T);

%Initialize inputs for afferent input A
onsetDelay  = pStruct.onsetDelay;        %Wait time to ripple start
% learnOL     = pStruct.learnOverlap;
learnOL     = round(pStruct.learnOverlap*learnDur);
rampLen1    = pStruct.rampLen1;           %Duration ramp length e.g. 1/5 or 1/2

AIn2(cueN,onsetDelay+1:onsetDelay+testDur)= Iexcit2;  %Test cue pulse
squarea             = Iexcit1*learnDur;     %Area under the square pulse

%Build Stimulation Protocols Based on RampTypeFlag
if learnDur ~= 0
for i = 1:N
    tmpStart                                    = onsetDelay+1+(i-1)*(learnDur)-(i-1)*(learnOL);
    tmpEnd                                      = tmpStart+learnDur;
    if pStruct.rampTypeFlag == 1                %FR IMA
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
    elseif pStruct.rampTypeFlag == 2            %DR IMA
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    elseif pStruct.rampTypeFlag == 6            %BR IMA
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    end
    if pStruct.rampTypeFlag == 3                %FR IP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);  %Calculate Single Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
    elseif pStruct.rampTypeFlag == 4            %DR IP
        maxI                                    = squarea/(learnDur - rampLen1); %Calculate Double Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    elseif pStruct.rampTypeFlag == 7            %BR IP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    end
end
end

%% Build voltage noise
if pStruct.noiseFlag == 0
    ANoise = zeros(N,T);
elseif pStruct.noiseFlag == 1
    ANoise = rand(N,T)*pStruct.noiseAmp - 0.5*pStruct.noiseAmp;
elseif pStruct.noiseFlag == 2
    ANoise = pinknoise(T,N)'*pStruct.noiseAmp;
end

% disp(['rampLen1: ', num2str(pStruct.rampLen1), '; learnDur: ', num2str(pStruct.learnDur), '; OL:', num2str(pStruct.learnOverlap)]);

%Learn Phase
for t=1:T
    for j = 1:N
        if (aLearn(j,t)   - tha(j)) > 0; aaLearn(j)   = aLearn(j,t)   - tha(j); else aaLearn(j)   = 0; end %Threshold activity under 'tha' to 0
        if (hLearn(j,t)   - thh(j)) > 0; hhLearn(j)   = hLearn(j,t)   - thh(j); else hhLearn(j)   = 0; end %Threshold activity under 'thh' to 0
        if (aLearn(j,t)   - thc(j)) > 0; ccLearn(j)   = aLearn(j,t)   - thc(j); else ccLearn(j)   = 0; end %Threshold activity under 'thc' to 0
    end
    
    if pStruct.simTypeFlag == 1
        %Standard Linear without Adapatation
        daLearn = AIn1(:,t) + ANoise(:,t) + achLearn.*(aaLearn*W(:,:,t))' - (hhLearn*H)' - eta1.*aLearn(:,t);
        dhLearn = (aaLearn*AH)' - HAuto.*hhLearn' - eta2.*hLearn(:,t);
    elseif pStruct.simTypeFlag == 2
        %Linear variant with Adaptation
        daLearn = AIn1(:,t) + ANoise(:,t) + achLearn.*(aaLearn*W(:,:,t))' - (hhLearn*H)' - eta1.*aLearn(:,t) + mu.*cLearn(:,t).*(Ek - aLearn(:,t));
        dhLearn = (aaLearn*AH)' - HAuto.*hhLearn' - eta2.*hLearn(:,t);
        dcLearn = gm.*ccLearn' - om.*cLearn(:,t);
    end

    aLearn(:,t+1) = aLearn(:,t)+daLearn*dt; %Update next time point of a
    hLearn(:,t+1) = hLearn(:,t)+dhLearn*dt; %Update next h time point
    
    if pStruct.simTypeFlag == 2
        %Adaptation update
        cLearn(:,t+1) = cLearn(:,t) + dcLearn;
    end

    %Weight update
    if pStruct.wtNoiseFlag == 0
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    elseif pStruct.wtNoiseFlag == 1
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    end

    W(:,:,t+1) = W(:,:,t)+dW*dt;   %Update next Weight time point
end
WLearn = W(:,:,end);    %Final Wt state after learn phase
% W(:,:,1) = WLearn;  %Reset first Wt Mat index for test phase learning

%Test Phase
for t=1:T2
    for j = 1:N
        if (aTest(j,t)- tha(j)) > 0; aaTest(j)= aTest(j,t)- tha(j); else aaTest(j)= 0; end
        if (hTest(j,t)- thh(j)) > 0; hhTest(j)= hTest(j,t)- thh(j); else hhTest(j)= 0; end
        if (aTest(j,t)- thc(j)) > 0; ccTest(j)= aTest(j,t)- thc(j); else ccTest(j)= 0; end %Threshold activity under 'thc' to 0
    end
    
    if pStruct.simTypeFlag == 1
        %Standard Linear without Adapatation
        daTest = AIn2(:,t) + ANoise(:,t) + achTest.*(aaTest*WLearn)' - (hhTest* H)' - eta1.*aTest(:,t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
        dhTest = (aaTest*AH)' - HAuto.*hhTest' - eta2.*hTest(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
    elseif pStruct.simTypeFlag == 2
        %Linear variant with Adaptation
        daTest = AIn2(:,t) + ANoise(:,t) + achTest.*(aaTest*WLearn)' - (hhTest*H)' - eta1.*aTest(:,t) + mu.*cTest(:,t).*(Ek - aTest(:,t));
        dhTest = (aaTest*AH)' - HAuto.*hhTest' - eta2.*hTest(:,t);
        dcTest = gm.*ccTest' - om.*cTest(:,t);
    end
    
    aTest(:,t+1) = aTest(:,t)+daTest*dt;
    hTest(:,t+1) = hTest(:,t)+dhTest*dt;
    
    if pStruct.simTypeFlag == 2
        %Adaptation update
        cTest(:,t+1) = cTest(:,t) + dcTest;
    end
        
%     %Weight update
%     if pStruct.wtNoiseFlag == 0
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
% %         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
%     elseif pStruct.wtNoiseFlag == 1
% %         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
%     end
% 
%     W(:,:,t+1) = W(:,:,t)+dW;   %Update next Weight time point
end
% WTest = W(:,:,end);     %Final Wt state after test phase

actCell     = {aLearn,aTest};
hactCell    = {hLearn,hTest};
inCell      = {AIn1,AIn2};
wOut        = {W_pre,WLearn};

end