function [fhandle,ax1] = plotCA1RipEx(ca3a,ca1a,ca1A,locsR,pStruct)
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
N = size(ca3a,1);

ca3Clr = autumn(N);
ca1Clr = cool(N+3);
if pStruct.rampTypeFlag == 'ctl'
    legStr = {'CA3 Pyr','Control CA1 Pyr'};
%     ca3Clr = cool(N);
elseif pStruct.rampTypeFlag == 0
    legStr = {'CA3 Pyr','Square CA1 Pyr'};
elseif pStruct.rampTypeFlag == 1
    legStr = {'CA3 Pyr','FR IMA CA1 Pyr'};
elseif pStruct.rampTypeFlag == 2
    legStr = {'CA3 Pyr','DR IMA CA1 Pyr'};
elseif pStruct.rampTypeFlag == 3
    legStr = {'CA3 Pyr','FR IP CA1 Pyr'};
elseif pStruct.rampTypeFlag == 4
    legStr = {'CA3 Pyr','DR IP CA1 Pyr'};
elseif pStruct.rampTypeFlag == 5
    legStr = {'CA3 Pyr','BR IMA CA1 Pyr'};
elseif pStruct.rampTypeFlag == 7
    legStr = {'CA3 Pyr','BR IP CA1 Pyr'};
end

%Plot activity
fhandle = figure;  hold on; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.38]);

ax1 = gca;
if pStruct.rampTypeFlag ~= 'ctl'
    rr = rectangle('Position',[pStruct.stimDelay,-27,pStruct.T-pStruct.stimDelay,54],'FaceColor',[0.9,0.9,1],'EdgeColor','none');
end

plot(ca3a(1,:),'Color',ca3Clr(1,:))
plot(ca1a(1,:),'Color',ca1Clr(1,:))
for i = 2:N
    plot(ca3a(i,:),'Color',ca3Clr(i,:))
%     if i == 1; plot(ARamp(pStruct.cueN,:).*50-26.5,'k'); end
%     plot(ARamp(i,:).*50-26.5,'Color',gc(i,:))
end
for i = 1:N
    plot(ca1a(i,:),'Color',ca1Clr(i,:))
end

plot(ca1A(1,:).*50-26.5,'k')

plot([0 1000],[10 10],'k--');
scatter(locsR,10*ones(1,length(locsR)),'k^','filled')
ylim([-30 30]); % ylabel('Activation'); 
xlim([0,1000]); % xlabel('Time (ms)');
legend(legStr,'Location','northeast')
set(gca,'FontSize',24,'fontname','times')

end