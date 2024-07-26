%%% Sample code for Attractor Model
%LKW 2/6/2022 also see Hasselmo et al., J. Neurosci. 1995.
%Compares simple adaptation current using simulated intracellular Calcium dynamics
%against non-adapting simulations with the same neuronal parameters and inputs

clearvars
close all

saveFlag = 0;
saveDir = 'Your\Dir\Here\';
sName = '_FR_IMA';

rampTypeFlag = 1;
W       = 0.0;        %Weight of a to a (recurrency); try 0.033
H       = 0.0;        %Weight of h to a (feedback inhibition); try 0.035
Wh      = 0.05;         %Weight of a to h (feedforward activation of IN)
Hh      = 0.003;        %Inhibition self feedback, try 0.003
tha     = 4;            %Threshold for a
thh     = 4;            %Thres1hold for h
thc     = 4;                %Threshold for voltage-dependent Ca-currents 
eta     = 0.01;         %Decay constant
T       = 1400;         %Time steps, must be even
Iexcit  = 1;          %Stim strength
inDur   = 1000;           %Stim duration
pDelay  = 200;
hiThresh= 10;           %Threshold of high firing likelihood
rampPerc= 0.25;
% Ionic Currents and related parameters 
Ek  = -10;              %Reversal Potential of Potassium
Ena = 70;               %Sodium
Ecl = 0;                %Chlorine
mu  = 0.01;             %Ca-dependent K-current
gm  = 0.001;            %Gamma; voltage-dependent Ca-currents
om  = 0.001;            %Omega; constant for diffusion of intracellular Ca

%Init Square vectors
AAdapt = zeros(1,T);   % Initialize input vector
aAdapt = zeros(1,T);
hAdapt = zeros(1,T);
cAdapt = zeros(1,T);

%Init Ramp vectors
ARamp = zeros(1,T);
aRamp = zeros(1,T);
hRamp = zeros(1,T);
cRamp = zeros(1,T);

%Initialize inputs
rampLen = round(inDur*rampPerc);
AAdapt(pDelay + 1:pDelay+inDur) = Iexcit;

if rampTypeFlag == 0
    ARamp(pDelay + 1:pDelay+inDur) = Iexcit;
    platPt = pDelay + inDur;
elseif rampTypeFlag == 1
    ARamp(pDelay + 1:pDelay + inDur)  = Iexcit;
    ARamp(pDelay + 1:pDelay + rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
    platPt = pDelay + inDur;
    figspecs = 'r-';
    legSpecs = {'FR IMA','Square','Adaptation \tau'};
elseif rampTypeFlag == 2
    ARamp(pDelay + 1:pDelay + inDur)  = Iexcit;
    ARamp(pDelay + 1:pDelay + rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
    ARamp(pDelay + inDur-rampLen+1:pDelay + inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
    platPt = pDelay + inDur - rampLen;
    figspecs = 'b-';
    legSpecs = {'DR IMA'};
elseif rampTypeFlag == 3
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(pDelay + 1:pDelay + inDur)  = IRamp;
    ARamp(pDelay + 1:pDelay + rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
    platPt = pDelay + inDur;
    figspecs = 'r--';
    legSpecs = {'FR IP'};
elseif rampTypeFlag == 4
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(inDur - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ARamp(pDelay + 1:pDelay + inDur)  = IRamp;
    ARamp(pDelay + 1:pDelay + rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
    ARamp(pDelay + inDur-rampLen+1:pDelay + inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
    platPt = pDelay + inDur - rampLen;
    figspecs = 'b--';
    legSpecs = {'DR IP'};
elseif rampTypeFlag == 5
    ARamp(pDelay + 1:pDelay + inDur)  = Iexcit;
    ARamp(pDelay + inDur-rampLen+1:pDelay + inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
    platPt = pDelay + inDur - rampLen;
    figspecs = 'c-';
    legSpecs = {'BR IMA'};
elseif rampTypeFlag == 6
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(pDelay + 1:pDelay + inDur)  = IRamp;
    ARamp(pDelay + inDur-rampLen+1:pDelay + inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
    platPt = pDelay + inDur - rampLen;
    figspecs = 'c--';
    legSpecs = {'BR IP'};
end

% Run Simulation

for t=1:T-1
    if (aAdapt(t) - tha) > 0; aaAdapt = aAdapt(t) - tha; else aaAdapt = 0; end
    if (hAdapt(t) - thh) > 0; hhAdapt = hAdapt(t) - thh; else hhAdapt = 0; end
    if (aAdapt(t) - thc) > 0; ccAdapt = aAdapt(t) - thc; else ccAdapt = 0; end
    if (aRamp(t)  - tha) > 0; aaRamp  = aRamp(t)  - tha; else aaRamp  = 0; end
    if (hRamp(t)  - thh) > 0; hhRamp  = hRamp(t)  - thh; else hhRamp  = 0; end
    if (aRamp(t)  - thc) > 0; ccRamp  = aRamp(t)  - thc; else ccRamp  = 0; end

    %Linear with Adaptation variant
    daAdapt = AAdapt(t) + W*aaAdapt - H*hhAdapt - eta*aAdapt(t) + mu.*cAdapt(t).*(Ek - aAdapt(t));
    dcAdapt = gm.*ccAdapt - om.*cAdapt(t);
    dhAdapt = Wh*aaAdapt - Hh*hhAdapt - eta*hAdapt(t);

    daRamp  = ARamp(t) + W*aaRamp - H*hhRamp - eta*aRamp(t) + mu.*cRamp(t).*(Ek - aRamp(t));
    dcRamp  = gm.*ccRamp - om.*cRamp(t);
    dhRamp  = Wh*aaRamp - Hh*hhRamp - eta*hRamp(t);

    %Neuronal activity update
    aAdapt(t+1) = aAdapt(t)+daAdapt;
    hAdapt(t+1) = hAdapt(t)+dhAdapt;
    
    aRamp(t+1)  = aRamp(t)+daRamp;
    hRamp(t+1)  = hRamp(t)+dhRamp;

    %Adaptation update (should remain 0 for cLin)
    cAdapt(t+1) = cAdapt(t) + dcAdapt;
    cRamp(t+1)  = cRamp(t) + dcRamp;
end

aMax = max(aAdapt);
aPlat = aAdapt(pDelay + inDur);
aSag = aMax - (aMax - aPlat)*2/3;  %Time constant is 2/3 of sag
aTC = find(aAdapt > aSag,1,'last') - find(aAdapt == aMax);
disp(['aTC = ', num2str(aTC)])

aRmpMax = max(aRamp);
aRmpPlat = aRamp(platPt);
aRampSag = aRmpMax - (aRmpMax - aRmpPlat)*2/3;  %Time constant is 2/3 of sag
aRampTC = find(aRamp > aRampSag,1,'last') - find(aRamp == aRmpMax);
disp(['aRampTC = ', num2str(aRampTC)])
%% Plotting optional
set(0,'DefaultLineLineWidth',2)

if max(aAdapt) > 0
    maxfig = max(aAdapt);
else
    maxfig = 50;
end
maxfig = inf;

%Main activity Fig
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.4, 0.3, 0.5]);
axis([1 T min(aAdapt)-10 maxfig]); axis square; hold on
figa=gcf; axa=gca; colormap(gray);
set(figa,'Name','a - activation','numbertitle','off');
xlabel('Time (ms)'); ylabel('Activation')
set(gca,'FontSize',20,'fontname','times')

axes(axa); 
% plot([0, T],[hiThresh,hiThresh],'k--','HandleVisibility','off')
plot(aRamp,figspecs);
plot(aAdapt,'k'); 
% plot([pDelay,pDelay+inDur],[-5 -5],'k','HandleVisibility','off')
plot(find(aAdapt > aSag,1,'last'),aSag,'ko')
plot(find(aRamp > aRampSag,1,'last'),aRampSag,'ko')
yyaxis right
plot(AAdapt,'-black')
plot(ARamp-0.1,figspecs);
ylim([-0.5, 7])
ylabel('Input Amplitude')
yticklabels([0 1])
legend(legSpecs)

% % Adaptation Ca Concentration figure
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.4, 0.28, 0.5]);
% hold on;
% plot(cAdapt, 'r')
% plot(cLin, 'b')
% legend('Ca Adapt','Ca No-Adapt','Fontsize',16)
% ylabel('Ca Adaptation Current');
% set(gca,'FontSize',20,'fontname','times')
% xlim([0 T])
% 
% % Adaptation Ca-Driven K-Current figure
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.55, 0.4, 0.28, 0.5]);
% hold on;
% plot(mu*cAdapt.*(Ek - aAdapt), 'r')
% plot(mu*cLin.*(Ek - cLin), 'b')
% legend('K Adapt','K No-Adapt','Location','Southeast','Fontsize',16)
% ylabel('Ca Driven K-Current');
% set(gca,'FontSize',20,'fontname','times')
% xlim([0 T])

% %Input figure
% figure; aae = gca;
% set(gcf,'Units','Normalized','OuterPosition', [0.25,0.25,0.28,0.2]);
% plot(AAdapt,'k-'); hold on;
% plot(ALin,'k--')
% legend('A Square', 'A Ramp','FontSize',16)
% xlim([0 5000]); ylim([0 max(ALin)+0.02])
% aae.XTick = [];
% ylabel('Input')
% set(gca,'FontSize',20,'fontname','times')

%% Model Comparisons
% Must first store or laod two aAdapt outputs from above simulation named
% as aMod1 and aMod2 (e.g. one simulation run with square (0% ramp) and one
% run with ramped input pulse

modCompsFlag = 0;
if modCompsFlag == 1
mod1Max = max(aMod1);
mod1Plat = aMod1(inDur+pDelay);
mod1Sag = mod1Max - (mod1Max - mod1Plat)*2/3;  %Time constant is 2/3 of sag
mod1TC = find(aMod1 > mod1Sag,1,'last') - find(aMod1 == mod1Max);
mod2Max = max(aMod2);
mod2Plat = aMod2(inDur+pDelay);
mod2Sag = mod2Max - (mod2Max - mod2Plat)*2/3;  %Time constant is 2/3 of sag
mod2TC = find(aMod2 > mod2Sag,1,'last') - find(aMod2 == mod2Max);

disp(['Model 1 TC: ',num2str(mod1TC)])
disp(['Model 2 TC: ',num2str(mod2TC)])

figure; plot(aMod1,'k')
hold on;
plot(aMod2,'r');
% plot(AAdapt+50,'k')
plot(find(aMod1 > mod1Sag,1,'last'),mod1Sag,'ko')
plot(find(aMod2 > mod2Sag,1,'last'),mod2Sag,'ko')
yyaxis right
plot(AMod1,'k')
plot(AMod2,'r')
ylim([-0.5, 7])
yticklabels([0 1])
legend('Model 1','Model 2')
set(gca,'FontSize',20,'fontname','times')
end
%% Saves
if saveFlag == 1
    saveas(figa,[saveDir,'activityPlot',sName],'png')
    saveas(figa,[saveDir,'activityPlot',sName],'fig')
    saveas(figa,[saveDir,'activityPlot',sName],'svg')
end