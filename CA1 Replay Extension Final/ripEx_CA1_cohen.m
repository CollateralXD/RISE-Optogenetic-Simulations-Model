function [dOut,meansOut,dOutIN,seqLen] = ripEx_CA1_cohen(actCell,hactCell,pStruct)
%Run Stats on pyramidals only for the CA3Net Ripple Extend code
%Input: actCell = cell containing the activity of pyramidals for ramp,
%square and control
%       hactCell = cell with activity of interneurons
% 
%Output: cohen's d of the comparison to control and raw mean and stdev of
%opto ramp
ca1PRmp = actCell{2}; ca1PCtl = actCell{3};
ca1Irmp = hactCell{2}; hControl = hactCell{3};
actThresh = pStruct.actThresh; 
pksR = []; locsR = []; pksC = []; locsC = [];
pksRIN = []; locsRIN = []; pksCIN = []; locsCIN = [];

for  i = 1:pStruct.N    %Find the first peak, or return no peak if flat line
    [pksTmp,locsTmp] = findpeaks(ca1PRmp(i,:));
    [pksINTmp,~] = findpeaks(ca1Irmp(i,:));
    if ~isempty(locsTmp) && max(pksTmp) > actThresh(i)    %If max activation above threshold
        pksR = [pksR max(pksTmp)];
        tmpInd = find(ca1PRmp(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsR = [locsR tmpInd];
        pksRIN = [pksRIN pksINTmp(1)]; 
        tmpInd = find(ca1Irmp(i,:) > actThresh(i),1);  %Get index of threshold crossing point
        locsRIN = [locsRIN tmpInd];
    end
    [pksTmp,locsTmp] = findpeaks(ca1PCtl(i,:));
    [pksINTmp,~] = findpeaks(hControl(i,:));
    if ~isempty(locsTmp) && ~isempty(pksINTmp) && max(pksTmp) > actThresh(i)
        pksC = [pksC max(pksTmp)];
        tmpInd = find(ca1PCtl(i,:) > actThresh(i),1);  %Get index of threshold crossing point
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
end

end