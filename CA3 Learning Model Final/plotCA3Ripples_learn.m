function [handle,ax1,ax2] = plotCA3Ripples_learn(actCell,hactCell,inCell,locsCell,pksCell,pStruct)
%%% 
%Inputs: 
% actCell 1x2 cell of pyramidal activation matrices in NxT format
% hactCell 1x2 cell of inhibitory activation matrices NxT format
% inCell 1x2 cell of inputs to the system as matrices NxT format
% cueN integer of the cued node for the simulation

%Outputs:
% handle = handle of the created figure
% axN    = handle to axis N of the figure
%%% 

set(0,'DefaultLineLineWidth',2)

aRamp = actCell{1}; aControl = actCell{2};
hRamp = hactCell{1};hControl = hactCell{2};
ARamp = inCell{1};  AControl = inCell{2};
locsR = locsCell{1};locsC = locsCell{2};
pksR = pksCell{1};  pksC  = pksCell{2};
N = pStruct.N;

%Plot activity
handle = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.95]);

subplot(2,1,1)
ax1 = gca;
% plot(hRamp(cueN,:),'b-.')
legCell(1) = {strcat(num2str(pStruct.rampPerc*100),'% Ramp Input')}; %'IN1'
ct = length(legCell)+1;
rc = cool(N);
sc = gray(N+2);
for i = 1:N
    plot(ARamp(i,:).*5,'Color',sc(i,:));  hold on;
    plot(aRamp(i,:),'Color',rc(i,:))
%     plot(hRamp(i,:),'b-.')
%     legCell(ct) = {strcat('R Pyr ',num2str(i))};  ct = ct+1;
%     legCell(ct) = {strcat('R IN ',num2str(i))};   ct = ct+1;
end
legCell(ct) = {strcat('Learn Pyr')};
% scatter(locsR,pksR*0,'k^','filled')
ylim([-inf,35]); ylabel('Activation'); 
xlim([0,pStruct.TLearn]); % xlabel('Time (ms)');
legend(legCell,'Location','northeast')
set(gca,'FontSize',24, 'XTickLabel',[],'fontname','times')

subplot(2,1,2)
ax2 = gca;
plot(AControl(pStruct.cueN,:).*5,'k'); hold on; 
% plot(hControl(cueN,:),'r-.')
legCell = {'Test Input'}; %'IN1'
ct = length(legCell)+1;
for i = 1:N
    plot(aControl(i,:),'Color',sc(i,:))
%     plot(hSquare(i,:),'r-.')
%     legCell(ct) = {strcat('S Pyr ',num2str(i))};  ct = ct+1;
%     legCell(ct) = {strcat('S IN ',num2str(i))};   ct = ct+1;
end
legCell(ct) = {strcat('Test Pyr')};
% scatter(locsS,pksS*0,'k^','filled')
xlim([0,1000]); xlabel('Time (ms)'); 
ylim([-inf,35]); ylabel('Activation'); 
legend(legCell,'Location','northeast')
set(gca,'FontSize',24,'fontname','times')

end