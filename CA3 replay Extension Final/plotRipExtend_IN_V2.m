function [fhandle,ax1] = plotRipExtend_IN_V2(aRamp,hRamp,ARamp,pStruct)
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
N = size(hRamp,1);

nc = jet(N);
rc = gray(N+5);
if pStruct.rampTypeFlag == 'ctl'
    rampStr = 'Control INs';
elseif pStruct.rampTypeFlag == 'sqr'
    rampStr = 'Square INs';
elseif pStruct.rampTypeFlag == 1
    rampStr = 'FR IMA INs';
elseif pStruct.rampTypeFlag == 2
    rampStr = 'DR IMA INs';
elseif pStruct.rampTypeFlag == 3
    rampStr = 'FR IP INs';
elseif pStruct.rampTypeFlag == 4
    rampStr = 'DR IP INs';
elseif pStruct.rampTypeFlag == 5
    rampStr = 'BR IMA INs';
elseif pStruct.rampTypeFlag == 7
    rampStr = 'BR IP INs';
end

%Plot activity
fhandle = figure;  hold on; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.05, 0.5, 0.38]);

ax1 = gca;
if pStruct.rampTypeFlag ~= 'ctl'
    rr = rectangle('Position',[pStruct.stimDelay,-27,pStruct.T-pStruct.stimDelay,100],'FaceColor',[0.9,0.9,1],'EdgeColor','none');
end

for i = 1:N
    plot(hRamp(i,:),'Color',nc(i,:))
    plot(aRamp(i,:),'Color',rc(i,:),'LineStyle',':')
    if i == 1; plot(ARamp(pStruct.cueN,:).*50-26.5,'k'); end
end
plot([0 1000],[10 10],'k--');
ylim([-27 35]); % ylabel('Activation'); 
xlim([0,1000]); % xlabel('Time (ms)');

if pStruct.rampTypeFlag ~= 'ctl'
    legend({rampStr,'Pyr'},'Location','northeast')
else
    legend({rampStr,'Pyr','Input'},'Location','northeast')
end
set(gca,'FontSize',24,'fontname','times')

end