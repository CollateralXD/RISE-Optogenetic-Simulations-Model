function [handle,ax1,ax2] = plotCA3Ripples_fast(actCell,hactCell,inCell,locsCell,pksCell,cueN)
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
N = size(aRamp,1);

%Plot activity
handle = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.95]);

subplot(2,1,1)
ax1 = gca;
plot(ARamp(cueN,:).*40+25,'k');  hold on; 
% plot(hRamp(cueN,:),'b-.')
legCell = {'Ramp Input'}; %'IN1'
ct = length(legCell)+1;
rc = cool(N);
for i = 1:N
    plot(aRamp(i,:),'Color',rc(i,:))
%     plot(hRamp(i,:),'b-.')
%     legCell(ct) = {strcat('R Pyr ',num2str(i))};  ct = ct+1;
%     legCell(ct) = {strcat('R IN ',num2str(i))};   ct = ct+1;
end
legCell(ct) = {strcat('Ramp Pyr')};
% scatter(locsR,pksR*0,'k^','filled')
ylim([-20,35]); ylabel('Activation'); 
xlim([0,1000]); % xlabel('Time (ms)');
legend(legCell,'Location','northeast')
set(gca,'FontSize',24, 'XTickLabel',[],'fontname','times')

subplot(2,1,2)
ax2 = gca;
plot(AControl(cueN,:).*40+25,'k'); hold on; 
% plot(hControl(cueN,:),'r-.')
legCell = {'Control Input'}; %'IN1'
ct = length(legCell)+1;
sc = gray(N+2);
for i = 1:N
    plot(aControl(i,:),'Color',sc(i,:))
%     plot(hSquare(i,:),'r-.')
%     legCell(ct) = {strcat('S Pyr ',num2str(i))};  ct = ct+1;
%     legCell(ct) = {strcat('S IN ',num2str(i))};   ct = ct+1;
end
legCell(ct) = {strcat('Control Pyr')};
% scatter(locsS,pksS*0,'k^','filled')
xlim([0,1000]); xlabel('Time (ms)'); 
ylim([-20,35]); ylabel('Activation'); 
legend(legCell,'Location','northeast')
set(gca,'FontSize',24,'fontname','times')

end