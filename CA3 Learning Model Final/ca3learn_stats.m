function [pks,locs,pkDelta,INDelta] = ca3learn_stats(pStruct,actTmp,hactTmp)
%%%
% Inputs: pStruct of parameters from CA3 learn code
%         actTmp: 1x2 cell array of Pyr activity from ca3Seq_learn.m
%         hactTmp: 1x2 cell array of IN activity from ca3Seq_learn.m
% Outputs: pks and locations of 
%%%
testAct = actTmp{2}; testHAct = hactTmp{2};
pks = []; locs = [];
pksIN = []; locsIN = [];
for i = 1:size(testAct,1)
    [pksTmp,locsTmp] = findpeaks(testAct(i,:));
    [pksINTmp,locsINTmp] = findpeaks(testHAct(i,:));
    %Add first nonzero, above threshold activation for each node to pksTest
    if ~isempty(locsTmp)
        if pksTmp(1) >= pStruct.actThresh(i)
            pks(i) = pksTmp(1); locs(i) = locsTmp(1); 
            pksIN(i) = pksINTmp(1); locsIN(i) = locsINTmp(1);
        end
    end
end
pkDelta = diff(locs);
INDelta = diff(locsIN);
end