function [fh1,fh2] = plotSumComp_learnfun(shufDat1,shufDat2,shufDat3,pStruct)
% Inputs: shuffled data structs containing durShuf, seqShuf and CIs
% pStruct with plotting parameters
% LKW 9/15/2021
nRamps = pStruct.nRamps;
yUp = pStruct.yUp;
yDn = pStruct.yDn;

if pStruct.WFC == 'FRC'
    lineSpecs = {'r-','r-o','r-d'};
elseif pStruct.WFC == 'DRC'
    lineSpecs = {'b-','b-o','b-d'};
elseif pStruct.WFC == 'BRC'
    lineSpecs = {'c-','c-o','c-d'};
elseif pStruct.WFC == 'FRP'
    lineSpecs = {'r--','r--o','r--d'};
elseif pStruct.WFC == 'DRP'
    lineSpecs = {'b--','b--o','b--d'};
elseif pStruct.WFC == 'BRP'
    lineSpecs = {'c--','c--o','c--d'};
end

patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;

fh1 = figure();
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aaa = gca;
axis square; hold on;
plot(shufDat1.meanRampD,lineSpecs{1});
plot(shufDat2.meanRampD,lineSpecs{2});
plot(shufDat3.meanRampD,lineSpecs{3});
patch(patchXs,[shufDat1.dnCId,fliplr(shufDat1.upCId)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[shufDat2.dnCId,fliplr(shufDat2.upCId)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[shufDat3.dnCId,fliplr(shufDat3.upCId)],'k','EdgeColor','none','FaceAlpha',0.1);
xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
if iscell(pStruct.legendCell)
    legend(pStruct.legendCell,'FontSize',16,'location',pStruct.legendLoc);
end
ylim([yDn(1) yUp(1)]); xlim([1 nRamps]); aaa.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
set(gca,'FontSize',24,'fontname','times')

fh2 = figure();
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aab = gca;
axis square; hold on;
plot(shufDat1.meanActN,lineSpecs{1});
plot(shufDat2.meanActN,lineSpecs{2});
plot(shufDat3.meanActN,lineSpecs{3});
patch(patchXs,[shufDat1.dnCIN,fliplr(shufDat1.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[shufDat2.dnCIN,fliplr(shufDat2.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1);
patch(patchXs,[shufDat3.dnCIN,fliplr(shufDat3.upCIN)],'k','EdgeColor','none','FaceAlpha',0.1);
xlabel('Ramp Percentage'); ylabel('Mean Sequence Length')
% if iscell(pStruct.legendCell)
%     legend(pStruct.legendCell,'FontSize',16,'location','east');
% end
ylim([yDn(2) yUp(2)]); xlim([1 nRamps]); aab.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
set(gca,'FontSize',24,'fontname','times')

end