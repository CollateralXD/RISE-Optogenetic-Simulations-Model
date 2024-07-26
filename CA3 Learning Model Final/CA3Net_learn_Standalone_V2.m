%% Sample Standalone CA3 Pattern Learning Code
% Simulates a learning and test phase for ripple activity in CA3 region
% under varying input waveform conditions
% LKW 8/16/21

clearvars
close all
% rng(5)      %For reproducibility

%1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 6 = BR IMA; 7 = BR IP
rampTypeFlag        = 1;
simTypeFlag         = 2;        %1 = Linear; 2 = Linear with Adaptation
noiseFlag           = 0;        %0 = no noise; 1 = White noise; 2 = 1/f or pink noise
wtNoiseFlag         = 0;        %0 = no noise; 1 = Lognormal noise distribution using parameters Mu and Sigma
saveFlag            = 0;
saveDir = 'Your\Dir\Here\';

disp(['Ramp Type ', num2str(rampTypeFlag),'; Sim Type ', num2str(simTypeFlag), '; Noise type ', num2str(noiseFlag), '; Wt Mat Noise type ', num2str(wtNoiseFlag)]);

cueN                = 1;            %Node index for cue
N                   = 15;           %# Nodes
Ww                  = 0.035;        %Baseline Pyr-Pyr weight strength try 0.035 for adaptation
Hh                  = 0.05;         %Baseline Weight strength IN-Pyr
Wh                  = 0.05;         %Baseline Weight strength Pyr-IN
HAuto               = 0.003;        %Baseline Weight IN-IN
tha                 = 4*ones(1,N);  %Activation Threshold Pyr
thh                 = 4*ones(1,N);  %Activation Threshold IN
thc                 = 4*ones(1,N);
eta1                = 0.01;         %Pyr Decay constant
eta2                = 0.01;         %IN Decay constant
noiseAmp            = 0.1*Ww;          %Center of lognormal weight noise if included.
noiseMu             = 1;            %Mean of lognormal wt distro
noiseSigma          = 0.4;            %Variance of lognorm wt distro
dt                  = 1;
T                   = 1500/dt;         %Time steps, must be even
Iexcit1             = 0.5;          %Learn phase afferent input
Iexcit2             = 1;            %Test phase afferent input
if wtNoiseFlag == 1
    achLearn        = 0.1;          %Ach state during learn phase; 0 = high Ach; 1 = low Ach
    e3              = 0.01;         %Learning rate constant
else
    achLearn        = 0.1;
    e3              = 0.001;
end
achTest             = 1;            %Ach state during test phase
e1                  = 0.002;         %Learning Decay Pre synaptic
onsetDelay          = 50/dt;           %Wait time for first pattern element to arise
learnDur            = 80/dt;           %Duration (ms) of learn pulses
OLPerc              = 0.60;           %Overlap (ms) of learn pattern elements
learnOL             = round(OLPerc*learnDur);
testDur             = 20/dt;           %Duration (ms) of test cue pulse
rampPerc            = 0.50;         %Percentage of ramp
rampLen1            = round(learnDur*rampPerc);  %Duration ramp length learn phase
AIn1                = zeros(N,T);   %Learn phase afferent input
AIn2                = zeros(N,T);   %Test phase
AIn2(cueN,onsetDelay:onsetDelay+testDur) = Iexcit2;  %Sq Pulse cue test phase
squarea             = Iexcit1*learnDur;     %Area under the square pulse

% Ionic Currents and related parameters 
Ek  = -10;              %Reversal Potential of Potassium
mu  = 0.01;             %Ca-dependent K-current
gm  = 0.001;            %Gamma; voltage-dependent Ca-currents
om  = 0.001;            %Omega; constant for diffusion of intracellular Ca

%Build Stimulation Protocols Based on RampTypeFlag
for i = 1:N
    tmpStart                                    = onsetDelay+1+(i-1)*(learnDur)-(i-1)*(learnOL);
    tmpEnd                                      = tmpStart+learnDur;
    if rampTypeFlag == 1                        %FRC
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
    elseif rampTypeFlag == 2                    %DRC
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    elseif rampTypeFlag == 6
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    end
    if rampTypeFlag == 3                        %FRP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);  %Calculate Single Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
    elseif rampTypeFlag == 4                    %DRP
        maxI                                    = squarea/(learnDur - rampLen1); %Calculate Double Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    elseif rampTypeFlag == 7                    %BRP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    end
end

%% Initialize Weight Matrices
if wtNoiseFlag == 0
    W  = zeros(N,N,T);    %Init wt mat
    AH = zeros(N,N);
    H  = zeros(N,N);
    
    for i = 1:N
        H(i,i)  = Hh;   %Direct feedback IN to Pyr
        AH(i,i) = Wh;   %Direct excitation pyr to IN
    end
    
elseif wtNoiseFlag == 1
    W = lognrnd(noiseMu,noiseSigma,N,N).*noiseAmp;          %Lognormal rand wt mat
    AH = zeros(N,N)*Wh/30;    %Wt mat of a to h (feedforward activation of IN)
    H  = zeros(N,N)*Hh/20;    %Weight of h to a (feedback inhibition)
    
    ahGain = 1;
    hGain = 1;
    for i = 1:N         %Build weight mats
        H(i,i)  = Hh;   %Direct feedback IN to Pyr
        AH(i,i) = Wh;   %Direct excitation pyr to IN
        if i <= N - 1   %Forward 1
            H(i,i+1)  = Hh/hGain;   %IN  2 Pyr
            AH(i,i+1) = Wh/ahGain;   %Pyr 2 IN
        end
        if i <= N - 2   %Forward 2
            H(i,i+2)  = Hh/hGain;   %IN  2 Pyr
            AH(i,i+2) = Wh/ahGain;   %Pyr 2 IN
        end
        if i <= N - 3   %Forward 3
            %         H(i,i+3)  = Hh/4;   %IN  2 Pyr
            %         AH(i,i+3) = Wh/4;   %Pyr 2 IN
        end
        if i > 1        %Back 1
            H(i,i-1)  = Hh/hGain;   %IN  2 Pyr
            AH(i,i-1) = Wh/ahGain;   %Pyr 2 IN
        end
        if i > 2        %Back 2
            H(i,i-2)  = Hh/hGain;   %IN  2 Pyr
            AH(i,i-2) = Wh/ahGain;   %Pyr 2 IN
        end
        if i > 3        %Back 3
            %         H(i,i-3)  = Hh/1;   %IN  2 Pyr
            %         AH(i,i-3) = Wh/4;   %Pyr 2 IN
        end
    end
end

% % Log Norm Rand Wt Distribution figure for manuscript (set N = 10000)
% W_sig0_5 = lognrnd(noiseMu,0.5,10000,1).*noiseAmp;
% W_sig1 = lognrnd(noiseMu,1,10000,1).*noiseAmp;
% W_sig0_25 = lognrnd(noiseMu,0.25,10000,1).*noiseAmp;
% figure; hold on;
% [~,edges] = histcounts(log10(W_sig1));
% h1 = histogram(W_sig0_25,10.^edges);
% h1.Normalization = 'probability';
% h2 = histogram(W_sig0_5,10.^edges);
% h2.Normalization = 'probability';
% h3 = histogram(W_sig1,10.^edges);
% h3.Normalization = 'probability';
% ylabel('Probability')
% xlabel('Weight Strength')
% legend('\sigma 0.25','\sigma 0.5','\sigma 1')
% set(gca, 'xscale','log','fontname','times','FontSize',20)

% Pre-learn wt figure for manuscript
figure; aac = gca; %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.3, 0.5]); hold on;
imagesc(flipud(W(:,:,1))); axis square; wcb = colorbar; %wcb.Label.String = "Weight Strength";
% xlabel("Post-Synaptic CA3"); ylabel("Pre-Synaptic CA3")
aac.YTick = 1:7:15; yticklabels(linspace(15,1,3));
aac.XTick = 1:7:15; xticklabels(linspace(1,15,3));
title("Learn Phase Start Network State")
set(gca,'fontname','times','FontSize',24)

%Initialize activity vectors for learn/test phases
aLearn = zeros(N,T);
hLearn = zeros(N,T);
cLearn = zeros(N,T);
aTest  = zeros(N,T);
hTest  = zeros(N,T);
cTest  = zeros(N,T);

aaL_e1 = zeros(N,N,T);
aaL_e2 = zeros(N,N,T);

%% Run learn block

for t=1:T
    %Neuron Activity
    for j = 1:N
        if (aLearn(j,t)   - tha(j)) > 0; aaLearn(j)   = aLearn(j,t)   - tha(j); else aaLearn(j)   = 0; end %Threshold activity under 'tha' to 0
        if (hLearn(j,t)   - thh(j)) > 0; hhLearn(j)   = hLearn(j,t)   - thh(j); else hhLearn(j)   = 0; end %Threshold activity under 'thh' to 0
        if (aLearn(j,t)   - thc(j)) > 0; ccLearn(j)   = aLearn(j,t)   - thc(j); else ccLearn(j)   = 0; end %Threshold activity under 'thc' to 0
    end
    
    if simTypeFlag == 1
        %Standard Linear without Adapatation
        daLearn = AIn1(:,t) + achLearn.*(aaLearn*W(:,:,t))' - (hhLearn*H)' - eta1.*aLearn(:,t);
        dhLearn = (aaLearn*AH)' - HAuto.*hhLearn' - eta2.*hLearn(:,t);
    elseif simTypeFlag == 2
        %Linear variant with Adaptation
        daLearn = AIn1(:,t) + achLearn.*(aaLearn*W(:,:,t))' - (hhLearn*H)' - eta1.*aLearn(:,t) + mu.*cLearn(:,t).*(Ek - aLearn(:,t));
        dhLearn = (aaLearn*AH)' - HAuto.*hhLearn' - eta2.*hLearn(:,t);
        dcLearn = gm.*ccLearn' - om.*cLearn(:,t);
    end
    
    aLearn(:,t+1) = aLearn(:,t)+daLearn*dt; %Update next time point of a
    hLearn(:,t+1) = hLearn(:,t)+dhLearn*dt; %Update next h time point
    
    if simTypeFlag == 2
        %Adaptation update
        cLearn(:,t+1) = cLearn(:,t) + dcLearn;
    end

    %Weight Update (suppressed by Ach state)
    if wtNoiseFlag == 0
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    elseif wtNoiseFlag == 1
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    end

    W(:,:,t+1) = W(:,:,t)+dW*dt;   %Update next Weight time point
end

%% Plot Learn
set(0,'DefaultLineLineWidth',2)
cPyr = autumn(N);
cIN  = jet(N);
cAIn = gray(N+2);

% % Plot learn afferent input
% figure; aaa = gca;
% imagesc(flipud(AIn1));
% % title('Afferent Input Pattern');
% acb = colorbar; 
% aaa.YTick = 1:7:15; yticklabels(linspace(15,1,3));
% set(gca,'fontname','times','FontSize',24)

% Plot learn activity
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.35]); hold on;
for i = 1:N
    plot(aLearn(i,:), 'Color', cPyr(i,:));
    plot(AIn1(i,:)*5, 'Color', cAIn(i,:));
end
xlabel('Time (ms)'); ylabel('Activation'); 
if rampTypeFlag == 1 && rampPerc == 0
    legend('Square','Input A Learn','location','southeast')
elseif rampTypeFlag == 1
    legend('FR, IMA','location','southeast')
elseif rampTypeFlag == 2
    legend('DR, IMA','location','southeast')
elseif rampTypeFlag == 3
    legend('FR, IP','location','southeast')
elseif rampTypeFlag == 4
    legend('DR, IP','location','southeast')
elseif rampTypeFlag == 6
    legend('BR, IMA','location','southeast')
elseif rampTypeFlag == 7
    legend('BR, IP','location','southeast')
end
% title('Learn Phase')
set(gca,'fontname','times','FontSize',24);

% Plot Post Learn Weights
figure; aaw = gca; %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.3, 0.5]); hold on;
imagesc(flipud(W(:,:,end))); axis square; wcb = colorbar; 
%     xlabel("Post-Synaptic CA3"); ylabel("Pre-Synaptic CA3")
%     wcb.Label.String = "Weight Strength";
aaw.YTick = 1:7:15; yticklabels(linspace(15,1,3));
aaw.XTick = 1:7:15; xticklabels(linspace(1,15,3));
title("Learn Phase Final Network State")
set(gca,'fontname','times','FontSize',24)

% Plot Post-learn - pre-learn
W_pre = W(:,:,1);
W_post = W(:,:,end);
figure; aaz = gca;
imagesc(flipud(W_post-W_pre)); axis square; wcb = colorbar;
aaz.YTick = 1:7:15; yticklabels(linspace(15,1,3));
aaz.XTick = 1:7:15; xticklabels(linspace(1,15,3));
title('Post - Pre Weights')
set(gca,'fontname','times','FontSize',24)

%% Run test block
W(:,:,1) = W(:,:,end);       %Set starting time point to learn mat for test phase

for t=1:T
    %Neuron Activity
    for j = 1:N
        if (aTest(j,t)   - tha(j)) > 0; aaTest(j)   = aTest(j,t)   - tha(j); else aaTest(j)   = 0; end %Threshold activity over 'tha' to 0
        if (hTest(j,t)   - thh(j)) > 0; hhTest(j)   = hTest(j,t)   - thh(j); else hhTest(j)   = 0; end %Threshold activity over 'thh' to 0
        if (aTest(j,t)   - thc(j)) > 0; ccTest(j)   = aTest(j,t)   - thc(j); else ccTest(j)   = 0; end %Threshold activity under 'thc' to 0
    end
    
    if simTypeFlag == 1
        %Standard Linear without Adapatation
        daTest = AIn2(:,t) + achTest.*(aaTest*W(:,:,t))' - (hhTest* H)' - eta1.*aTest(:,t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
        dhTest = (aaTest*AH)' - HAuto.*hhTest' - eta2.*hTest(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
    elseif simTypeFlag == 2
        %Linear variant with Adaptation
        daTest = AIn2(:,t) + achTest.*(aaTest*W(:,:,t))' - (hhTest*H)' - eta1.*aTest(:,t) + mu.*cTest(:,t).*(Ek - aTest(:,t));
        dhTest = (aaTest*AH)' - HAuto.*hhTest' - eta2.*hTest(:,t);
        dcTest = gm.*ccTest' - om.*cTest(:,t);
    end

    aTest(:,t+1) = aTest(:,t)+daTest*dt; %Update next time point of a
    hTest(:,t+1) = hTest(:,t)+dhTest*dt; %Update next h time point
    
    if simTypeFlag == 2
        %Adaptation update
        cTest(:,t+1) = cTest(:,t) + dcTest;
    end

    %Weight Update (suppressed by Ach state)
    if wtNoiseFlag == 0
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    elseif wtNoiseFlag == 1
%         dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Standard Hebbian with Ceiling at Ww
        dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn') - (e1*W(:,:,t).*repmat(aaLearn,[N,1])); %LearnRate * hebbian - presynaptic-gated decay
    end
        
    W(:,:,t+1) = W(:,:,t)+dW*dt;   %Update next Weight time point
end

%% Calculate Peaks and Locs of threshold crossing
pksT = zeros(N,1); locsT = zeros(N,1);
actThresh = 10*ones(N,1);
for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(aTest(i,:));
    if ~isempty(locsTmp); pksT(i) = pksTmp(1); locsT(i) = locsTmp(1); end
end
ttt = [];
for j = 1:numel(pksT)     %For each peak in the sim
    indtmp = find(aTest(j,:)>actThresh(j),1);
    if ~isempty(indtmp)
        ttt = [ttt,indtmp];
    end
end
    
%% Plot Test options

% % Test Phase Afferent input figure for manuscript
% figure; aaa = gca;
% imagesc(flipud(AIn2));
% xlabel('Time (ms)'); ylabel('Node #'); title('Afferent Input Pattern');
% aaa.YTick = 1:7:15; yticklabels(linspace(15,1,3));
% acb = colorbar; acb.Label.String = "Afferent Strength";
% set(gca,'fontname','times','FontSize',20)
% breakxaxis([250 1250],0.03);

aat = figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.35]);
% subplot(2,1,1)
for i = 1:N
    plot(aTest(i,:),'Color',cPyr(i,:)); hold on;
    plot(AIn2(i,:)*5, 'Color', cAIn(i,:));
end
plot([0,T-1],[10,10],'k--','LineWidth',1)
scatter(ttt,ttt*0+10,'^','k','filled');
% xlabel('Time (ms)'); ylabel('Activation'); 
if rampTypeFlag == 1 && rampPerc == 0
    legend('Square Pyramidal','Test Cue','location','southeast')
elseif rampTypeFlag == 1
    legend('FR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 2
    legend('DR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 3
    legend('FR IP Pyramidal','location','southeast')
elseif rampTypeFlag == 4
    legend('DR IP Pyramidal','location','southeast')
elseif rampTypeFlag == 6
    legend('BR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 7
    legend('BR IP Pyramidal','location','southeast')
end
% title('Test Phase Pyrs')
xlim([0,T-1]); ylim([-30,25])
set(gca,'fontname','times','FontSize',24);

%% Saves
if saveFlag == 1
    
disp('saving vars and figs')
OLStr = ['LDur_',num2str(learnDur),'_OL_',num2str(100*OLPerc),'_RmpPct_',num2str(100*rampPerc)];

if rampTypeFlag == 1 && rampLen1 == 0
    sBase = ['Square_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 1
    sBase = ['FR_IMA_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 2
    sBase = ['DR_IMA_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 6
    sBase = ['BR_IMA_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 3
    sBase = ['FR_IP_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 4
    sBase = ['DR_IP_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
elseif rampTypeFlag == 7
    sBase = ['BR_IP_',OLStr];
    saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
    saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
    saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
end
end
