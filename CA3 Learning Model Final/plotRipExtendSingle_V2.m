function [fhandle,ax1] = plotRipExtendSingle_V2(aRamp,ARamp,locsR,pStruct)
%%% Sample code for Replay Extension Model
%LKW 8/19/21
%%% 
%Inputs: 
% aRamp = matrix of pyramidal activations in NxT format
% hRamp = matrix of interneuron activities in NxT format
% ARamp = matrix of afferent pulse values in NxT format
% cueN  = double integer of the cued node for the simulation

%Outputs:
% handle = handle of the created figure
% axN    = handle to axis N of the figure
%%% 

set(0,'DefaultLineLineWidth',2)
N = size(aRamp,1);

rc = autumn(N);
if pStruct.rampTypeFlag == 'ctl'
    rampStr = 'Control Pyramidal';
    rc = cool(N);
elseif pStruct.rampTypeFlag == 'sqr'
    rampStr = 'Square Pyramidal';
elseif pStruct.rampTypeFlag == 1
    rampStr = 'FR IMA Pyramidal';
elseif pStruct.rampTypeFlag == 2
    rampStr = 'DR IMA Pyramidal';
elseif pStruct.rampTypeFlag == 3
    rampStr = 'FR IP Pyramidal';
elseif pStruct.rampTypeFlag == 4
    rampStr = 'DR IP Pyramidal';
elseif pStruct.rampTypeFlag == 5
    rampStr = 'BR IMA Pyramidal';
elseif pStruct.rampTypeFlag == 7
    rampStr = 'BR IP Pyramidal';
end

%Plot activity
fhandle = figure;  hold on; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.38]);

ax1 = gca;
if pStruct.rampTypeFlag ~= 'ctl'
    rr = rectangle('Position',[pStruct.stimDelay,-27,pStruct.T-pStruct.stimDelay,54],'FaceColor',[0.9,0.9,1],'EdgeColor','none');
end

for i = 1:N
    plot(aRamp(i,:),'Color',rc(i,:))
    if i == 1; plot(ARamp(pStruct.cueN,:).*50-26.5,'k'); end
end
plot([0 1000],[10 10],'k--');
scatter(locsR,10*ones(1,length(locsR)),'k^','filled')
ylim([-27 27]); % ylabel('Activation'); 
xlim([0,1000]); % xlabel('Time (ms)');

if pStruct.rampTypeFlag ~= 'ctl'
    legend(rampStr,'Location','northeast')
else
    legend(rampStr,'Input','Location','northeast')
end
set(gca,'FontSize',24,'fontname','times')

end