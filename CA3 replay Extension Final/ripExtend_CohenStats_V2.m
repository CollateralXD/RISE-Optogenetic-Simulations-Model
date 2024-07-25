function [dOut,meansOut,dOutIN,seqLen] = ripExtend_CohenStats_V2(actCell,hactCell,pStruct)
%Run Stats on pyramidals only for the CA3Net Ripple Extend code
%Input: actCell = cell containing the activity of pyramidals for ramp,
%square and control
%       hactCell = cell with activity of interneurons
% 
%Output: cohen's d of the comparison to control and raw mean and stdev of
%opto ramp
aRamp = actCell{1}; aControl = actCell{2};
hRamp = hactCell{1}; hControl = hactCell{2};
actThresh = pStruct.actThresh; 
pksR = []; locsR = []; pksC = []; locsC = [];
pksRIN = []; locsRIN = []; pksCIN = []; locsCIN = [];

for  i = 1:pStruct.N    %Find the first peak, or return no peak if flat line
    [pksTmp,locsTmp] = findpeaks(aRamp(i,:));
    [pksINTmp,~] = findpeaks(hRamp(i,:));
    if ~isempty(locsTmp) && max(pksTmp) > actThresh(i)    %If max activation above threshold
        pksR = [pksR max(pksTmp)];
        tmpInd = find(aRamp(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsR = [locsR tmpInd];
        pksRIN = [pksRIN pksINTmp(1)]; 
        tmpInd = find(hRamp(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsRIN = [locsRIN tmpInd];
    end
    [pksTmp,locsTmp] = findpeaks(aControl(i,:));
    [pksINTmp,~] = findpeaks(hControl(i,:));
    if ~isempty(locsTmp) && ~isempty(pksINTmp) && max(pksTmp) > actThresh(i)
        pksC = [pksC max(pksTmp)];
        tmpInd = find(aControl(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsC = [locsC tmpInd];
        pksCIN = [pksCIN pksINTmp(1)]; 
        tmpInd = find(hControl(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsCIN = [locsCIN tmpInd];
    end
end

seqLen = numel(locsR);  % Get sequence length of opto extension before truncating
postOptoLocs = locsR > pStruct.stimDelay;  %Only use the peaks after opto stim onset
postOptoLocsIN = locsRIN > pStruct.stimDelay;
locsR = locsR(postOptoLocs);
locsRIN = locsRIN(postOptoLocsIN);

%Prep for stats
rampDeltas = diff(locsR);
rampINDeltas = diff(locsRIN);
ctlDeltas = diff(locsC);
ctlINDeltas = diff(locsCIN);

meansOut = mean(rampDeltas);
% sdsOut = std(rampDeltas);

try
    dOut = computeCohen_d(rampDeltas,ctlDeltas,'independent');
    dOutIN = computeCohen_d(rampINDeltas,ctlINDeltas,'independent');
catch 
    dOut = NaN;
    dOutIN = NaN;
end

end